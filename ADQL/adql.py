import math
import numpy as np
import sqlparse
import re

import astropy.units as u
from astropy.coordinates import SkyCoord

from spatial_index import SpatialIndex


class ADQL:

    """
    ADQL class

    ADQL is a variant of SQL that is understood by no DBMSs, the main
    difference being support for sky spatial constraints.

    One approach to supporting ADQL is to extend the DBMS with these
    ADQL-specific constructs.  The alternate approach taken here is
    to translate the ADQL into local SQL before sending it to the
    DBMS.

    Since we also want to make use of spatial indexing to speed up
    spatial searches, we include in this translation additional
    spatial index column constraints.  This requires reference to
    special added columns.

    The hope is that this code can be leveraged by others who either
    wish to add spatial indexing to their DBMS the same way we have
    or would like a starting place for developing code around their
    own DBMS.

    This package is complementary to but independent from our
    implemetation of the IVOA Table Access Protocol (TAP); a web
    service for managing remote database requests for astronomy
    in a general way.  TAP manages requests/results and is generally
    fed ADQL.  This ADQL translator is used by our TAP implementation
    to turn ADQL into SQL the local DBMS can handle directly.
    """

    def __init__(self, mode=SpatialIndex.HTM, level=7, encoding=None,
                 debugfile=None, racol='ra', deccol='dec',
                 xcol='x', ycol='y', zcol='z', indxcol=None):

        """
        ADQL() initialization allows you to set a number of parameters
        associated with the indexing, all of which have defaults.

        Parameters
        ----------
        mode : integer, optional, default=SpatialIndex.HTM
            The spatial indexing supports both Heirarchical Triangular Mesh
            (HTM) and Hierarchical Equal Area isoLatitude Pixelization
           (HEALPix) tesselation of the sphere.
        level : integer, optional, default=7
            Both tesselation schemes involve a z-ordered nesting of cells.
            This parameter defines the depth used for this database table.
        xcol : string, optional, default='x'
            The algorithm requires four special columns in the table: the
            sky position (RA, Dec) three-vector (x,y,z) coordinates (doubles)
            and the spatial index (integer).  This parameter is the column
            name for the sky position three-vector X component.
        ycol : string, optional, default='y'
            Column name for sky position three-vector Y component.
        zcol : string, optional, default='z'
            Column name for sky position three-vector Z component.
        indxcol : string, optional, default depends on mode and encoding
            though we recommend using a name based on mode and level, such
            as 'htm20' or 'hpx14'.
            Column name for spatial index.
        encoding : integer, optional, default depends on level
            Since the tesselations are a nesting of four values, one
            encoding is essentially "base 4" (e.g. 13112032) though the number
            is still stored as a decimal integer.  Early on twenty years ago
            we stored the indices this way as it made debugging easier and we
            still support this encoding.  We would not recommend it for any
            new database. Use SpatialIndex.BASE10 instead.
        debugfile : string, optional, default None
            If set, defines the file where debug messages will go.
        """

        self.debug = False

        self.depth = ''

        self.closing = False

        self.adql_tokens = []

        self.funcData  = []
        self.funcLevel = 0

        self.mode     = mode
        self.level    = level
        self.xcol     = xcol
        self.ycol     = ycol
        self.zcol     = zcol

        self.indxcol  = indxcol
        self.encoding = encoding

        if indxcol is None:

            if encoding is None:
                self.encoding = SpatialIndex.BASE4
                self.indxcol = 'spt_ind'

            else:
                if(mode == SpatialIndex.HTM):
                    self.indxcol = 'htm' + str(level)
                else:
                    self.indxcol = 'hpx' + str(level)

        if self.encoding is None:
            self.encoding = SpatialIndex.BASE10

        if debugfile is not None and debugfile != '':
            self.debug = True
            self.debugfile = open(debugfile, 'w+')

        self.fconvert = {'icrs'          : 'icrs',
                         'fk5'           : 'fk5',
                         'fk4'           : 'fk4',
                         'j2000'         : 'fk5',
                         'b1950'         : 'fk4',
                         'ecliptic'      : 'geocentrictrueecliptic',
                         'galactic'      : 'galactic',
                         'supergalactic' : 'supergalactic'}


    def _parseADQL(self, stmt):

        # This function doesn't actually parse the AQDL; that has
        # already been done.  Instead it walks through the parse
        # tree, collecting information needed to translate the
        # bits of AQDL we want to change into different local SQL.
        # _parsADQL() calls itself recursively when it encounters
        # tokens that are actually branch nodes in the ADQL tree.
        #
        # Since our output is still SQL, it is pretty much the
        # same as the input.  So here we build a sequential list
        # of the token strings, replacing the collections that
        # change with a placeholder ('GEOM') and populating
        # arrays of functions and arguments we will use to build
        # the new output when we put the SQL back together.
        #
        # A lot of this is assocated with the geometric functions;
        # We want to processes 'AND CONTAINS()' blocks specially,
        # so if we see an Identifier with that token value, we
        # trigger parameter collection.  While this special
        # processing mode is on, we look for Identifiers
        # 'POINT', 'CIRCLE', and 'POLYGON' and set flags
        # for these, too.  And finally, if in one of these
        # we see an IdentifierList, we capture the list as
        # the arguments to the function:
        #
        #    CONTAINS -|- POINT   -|- arg1
        #              |           |- arg2
        #              |
        #              |- CIRCLE  -|- arg1
        #              |           |- arg2
        #              |           |- arg3
        #              |
        #              |- POLYGON -|- args
        #                          |- arg2
        #                          |- ..
        #                          |- argn
        #
        # We do the same thing for DISTANCE().
        #
        # We also need to watch for 'TOP <n>' constructs since
        # ADQL (and some DBMSs) use this to limit the number
        # of records to return but other DBMSs use other
        # constructs (and not in the same place in the query).
        # For our purpose (Oracle, this constraint goes in the
        # WHERE clause as an additional 'AND ROWNUM <= <n>'
        # so we will have to figure out where to insert it.
        #
        # ADQL (and some DBMSs) use this to limit the number
        # of records to return but other DBMSs use other
        # constructs (and not in the same place in the query).
        # For our purpose (Oracle) this constraint goes in the
        # WHERE clause as an additional 'AND ROWNUM <= <n>'
        # so we will have to figure out where to insert it.
        #
        # However, because of the way ROWNUM works, there is
        # no simple way to implement the ADQL 'OFFSET' directive.

        old_depth = self.depth

        self.depth = self.depth + '   '

        for i in range(len(stmt.tokens)):

            token = stmt.tokens[i]

            # IDENTIFIER

            # Almost all of our analysis for special function translation is
            # done by looking a sequences of TOKENs below.  However the
            # Identifier is the best place to handle "table.column" constructs.

            if(isinstance(token, sqlparse.sql.Identifier)):
                if(self.debug):
                    self.debugfile.write(self.depth + 'Identifier: [' + str(token) + ']\n')

                if self.funcLevel == 2 and len(token.tokens) == 3 and str(token.tokens[1]) == '.':
                    self.args.append('DOT')

            if(self.debug):


                # STATEMENT

                # By the time we get here we should already be inside a
                # Statement but we'll keep this here for debugging purposes.

                if(isinstance(token, sqlparse.sql.Statement)):
                    self.debugfile.write(self.depth + 'Statement:\n')


                # COMMENT

                # Our use case doesn't allow for comments but we'll leave this
                # in for debugging.

                elif(isinstance(token, sqlparse.sql.Comment)):
                    self.debugfile.write(self.depth + 'Comment:\n')


                # WHERE

                # The where clause is often a combination of "contains()"
                # and standard SQL constrains but the "where" itself
                # should just be passed along.

                elif(isinstance(token, sqlparse.sql.Where)):
                    self.debugfile.write(self.depth + 'Where:\n')


                # IDENTIFIERLIST

                # If we are in the middle of processing a "contains()"
                # function the IdentifierList is probably the function
                # argument list.

                elif(isinstance(token, sqlparse.sql.IdentifierList)):
                    self.debugfile.write(self.depth + 'IdentifierList:\n')


                # TOKENLIST

                # Like PARENTHESIS, we can let the TOKEN logic deal with
                # things that shouldn't be output because of "contains()".

                elif(isinstance(token, sqlparse.sql.TokenList)):
                    self.debugfile.write(self.depth + 'TokenList:\n')


                # IDENTIFIER

                # The identifiers we want for our special processing are also
                # single TOKENs as well.

                elif(isinstance(token, sqlparse.sql.Identifier)):
                    self.debugfile.write(self.depth + 'Identifier:\n')


                # COMPARISON

                # There shouldn't be any comparisons inside the "contains()"
                # sections (though one right after them) so we can just
                # process these normally.

                elif(isinstance(token, sqlparse.sql.Comparison)):
                    self.debugfile.write(self.depth + 'Comparison:\n')


                # PARENTHESIS

                # We won't want the parentheses inside the "contains()"
                # processing to be reproduced in the output query but like
                # a lot of these we need to keep the tokens otherwise.
                # We'll handle this by disallowing the saving of tokens
                # in the generic TOKEN processing below under the right
                # conditions.

                elif(isinstance(token, sqlparse.sql.Parenthesis)):
                    self.debugfile.write(self.depth + 'Parenthesis:\n')

            if(    isinstance(token, sqlparse.sql.Statement)
                or isinstance(token, sqlparse.sql.Comment)
                or isinstance(token, sqlparse.sql.Where)
                or isinstance(token, sqlparse.sql.IdentifierList)
                or isinstance(token, sqlparse.sql.Identifier)
                or isinstance(token, sqlparse.sql.TokenList)
                or isinstance(token, sqlparse.sql.Comparison)
                or isinstance(token, sqlparse.sql.Parenthesis)):
                    self._parseADQL(token)


            # TOKEN

            # Almost all of the function parameter collection is here.

            # If we are handling "contains()" and encounter a token,
            # it should either be a closing parenthesis (ending some
            # part of our functions) or a token we don't need (like
            # whitespace or commas).

            elif(not isinstance(token, sqlparse.sql.Token)):

                if(self.debug):
                    self.debugfile.write(self.depth + 'ERROR: Not a token: [' + str(token) + ']\n')

            elif(isinstance(token, sqlparse.sql.Token)):

                if(self.debug):
                    self.debugfile.write(self.depth + '[' + str(token) + ']\n')


                # Check through all the TOKENS.  Look for stuff associated with
                # the special ADQL "geometry" functions.

                # First check for Closing parenthesis which triggers
                # level change

                try:
                    if(token.value == ')'):

                        # Right parenthesis closing current 'args' array.
                        # Time to copy it to func['args'].

                        if(self.funcLevel == 2):
                            self.func['args'] = self.args
                            self.gargs.append(self.func)
                            self.funcLevel = 1

                        # Right parenthesis closing current 'args' array.
                        # Time to copy it to func['args'].

                        elif(self.funcLevel == 1):
                            self.gfunc['args'] = self.gargs
                            self.funcData.append(self.gfunc)
                            self.funcLevel = 0
                            self.closing = True

                except Exception:
                    pass


                # Level 0: check for primary functions

                if(self.funcLevel == 0):

                    if(token.value == 'distance'):
                        self.gfunc = {}
                        self.gfunc['name'] = 'distance'
                        self.gargs = []
                        self.funcLevel = 1
                        self.adql_tokens.append('GEOM')

                    elif(token.value == 'contains'):
                        self.gfunc = {}
                        self.gfunc['name'] = 'contains'
                        self.gargs = []
                        self.funcLevel = 1
                        self.adql_tokens.append('GEOM')

                    elif(token.value == 'jcg'):
                        self.gfunc = {}
                        self.gfunc['name'] = 'jcg'
                        self.gargs = []
                        self.funcLevel = 1
                        self.adql_tokens.append('GEOM')

                    elif(token.value != ')'):
                        self.adql_tokens.append(str(token))
                        self.closing = False

                    elif(token.value == ')' and self.closing is False):
                        self.adql_tokens.append(str(token))
                        self.closing = False


                # Level 1: check for support functions

                elif(self.funcLevel == 1):

                    if(token.value == 'point'):
                        self.func = {}
                        self.func['name'] = 'point'
                        self.args = []
                        self.funcLevel = 2

                    elif(token.value == 'circle'):
                        self.func = {}
                        self.func['name'] = 'circle'
                        self.args = []
                        self.funcLevel = 2

                    elif(token.value == 'box'):
                        self.func = {}
                        self.func['name'] = 'box'
                        self.args = []
                        self.funcLevel = 2

                    elif(token.value == 'polygon'):
                        self.func = {}
                        self.func['name'] = 'polygon'
                        self.args = []
                        self.funcLevel = 2

                    elif(token.value == 'jcg'):
                        self.func = {}
                        self.func['name'] = 'jcg'
                        self.gargs.append(token.value)


                # Level 2: Finally, all non-punctuation/whitespace
                # are added as level-2 arguments.

                elif(self.funcLevel == 2):

                    if(token.value != '(' \
                        and str(token.ttype)       != 'Token.Punctuation' \
                        and str(token.ttype)[0:21] != 'Token.Text.Whitespace' \
                        and str(token.ttype)[0:14] != 'Token.Operator' \
                        and str(token.ttype)       != 'Token.Comparison' \
                        and str(token.ttype)[0:13] != 'Token.Keyword'):
                            self.args.append(token.value)

        self.depth = old_depth


    def _geomConstraint(self, i):

        # This is where we convert the "contains()" function argument
        # into constraints our DBMS can use directly.  There are
        # multiple possibilities, some of which we cannot handle at
        # all and for which we will error off, some where we can
        # construct a constrain but which cannot make use of our
        # spatial index, and so where we can include spatial index
        # constraints and which will therefore run very fast.
        #
        # Basically, "contains()" requires that one geometric shape
        # be completely contained in another.  The shapes can be points,
        # circles (cones on the sky), or polygons (with 'box' as a
        # special case).  These shapes can be defined using fixed
        # values or database column names.
        #
        # Below we present a number of examples and explain what we
        # can do with each one.  For some of these, we use 'circle'
        # as a shorthand for circle or polygon or box.
        #
        #
        #    contains(point('icrs', ra, dec), circle('icrs', 45., 32., 0.5))
        #
        #    A point inside a fixed circle or convex polygon is the
        #    most important example for a couple of reasons.  First,
        #    it is by far the most common scientifically and second,
        #    it is the one only one where we can make effective use
        #    of our spatial index.
        #
        #
        #    contains(circle('icrs', ra, dec, 0.3), circle('icrs', 45., 32., 0.5))
        #
        #    A circle with a fixed radius inside another fixed region
        #    can also be indexed since that first fixed radius can
        #    be transferred to the second region mathematically.
        #    This may seem like it would not be of much use but in
        #    fact it lets us approximate extended objects (like small
        #    images) by their center and bounding circle.  For instance,
        #    all 2MASS images have a diagonal size of 0.134 degrees.
        #    A query of this sort will return a few extra objects but
        #    these are easily remove by post-filtering.
        #
        #
        #    contains(circle('icrs', ra, dec, radius), circle('icrs', 45., 32., 0.5))
        #    contains(point('icrs', 45., 32.), circle('icrs', ra, dec, radius))
        #
        #    However, when the radius is variable (and especially
        #    when there is no useful maximum) we can construct a
        #    constraint as a filter (see below) but cannot include
        #    any spatial index constraints.
        #
        #
        #    contains(polygon('icrs', ra1, dec1, ra2, dec2, ra3, dec3, ra4, dec4),
        #             circle('icrs', 45., 30., 0.5))
        #
        #    Anything beyond the above may make sense arithmetically
        #    but cannot be cast into simple SQL.  If there is a
        #    enough of a demand for it, we would need to extend the
        #    DBMS engine itself, something we want to avoid as part
        #    of this package.

        funcName = self.funcData[i]['name']


        # CONTAINS FUNCTION

        if funcName == 'contains':

            func1 = self.funcData[i]['args'][0]['name']
            func2 = self.funcData[i]['args'][1]['name']

            args1 = self.funcData[i]['args'][0]['args']
            args2 = self.funcData[i]['args'][1]['args']

            val   = self.funcData[i]['val']


            # Point in circle

            spt = SpatialIndex()

            nargs1 = len(args1)
            nargs2 = len(args2)

            geomstr = ''

            if( func1 == 'point'
            and nargs1 == 3
            and self._is_number(args1[1]) is False
            and self._is_number(args1[2]) is False

            and func2 == 'circle'
            and nargs2 == 4
            and self._is_number(args2[1]) is True
            and self._is_number(args2[2]) is True
            and self._is_number(args2[3]) is True):

                coordsys = ''
                try:
                    coordsys = self._frame_lookup(args1[0])
                except Exception as e:
                    raise Exception(str(e) + ' in point() function.\n')

                racol  = args1[1]
                deccol = args1[2]

                if(racol != 'ra'):
                    raise Exception('Invalid ra column reference in point() function.\n')

                if(deccol != 'dec'):
                    raise Exception('Invalid dec column reference in point() function.\n')


                try:
                    coordsys = self._frame_lookup(args2[0])
                except Exception as e:
                    raise Exception(str(e) + ' in circle function.\n')

                lon = float(args2[1])
                lat = float(args2[2])

                if(lat < -90.):
                    raise Exception('Invalid latitude value in circle() function.')

                if(lat > 90.):
                    raise Exception('Invalid latitude value in circle() function.')

                radius = float(args2[3])

                if(radius < 0.):
                    raise Exception('Invalid radius value in circle() function.')

                if(radius > 180.):
                    raise Exception('Invalid radius value in circle() function.')

                try:
                    coord = SkyCoord(lon * u.degree, lat * u.degree, frame=coordsys)
                except Exception as e:
                    raise Exception(str(e))

                ra  = coord.icrs.ra.deg
                dec = coord.icrs.dec.deg

                try:
                    retval = spt.cone_search(ra, dec, radius, self.mode, \
                        self.level, self.xcol, self.ycol, self.zcol, \
                        self.indxcol, self.encoding)

                except Exception as e:
                    raise Exception(str(e))

                if(val == '0'):
                    retval['geom_constraint'] = retval['geom_constraint'].replace('<=', '>')
                    geomstr = "(" + retval['geom_constraint'] + ")"

                else:
                    geomstr = "(" + retval['geom_constraint'] + ") AND (" + retval['index_constraint'] + ")"

            # Point in polygon

            elif(func1 == 'point'
             and nargs1 == 3
             and self._is_number(args1[1]) is False
             and self._is_number(args1[2]) is False

    #         and func2 == 'polygon' and nargs2 >= 7):
             and func2 == 'polygon'):

                if (nargs2 < 7):
                    raise Exception('Polygon specification requires at least 6 numbers.')

                try:
                    coordsys = self._frame_lookup(args1[0])
                except Exception as e:
                    raise Exception(str(e) + ' in point() function.')

                try:
                    coordsys = self._frame_lookup(args2[0])
                except Exception as e:
                    raise Exception(str(e) + ' in polygon() function.')

                if(nargs2 % 2 == 0):
                    raise Exception('Polygon specification requires an even number of values (>=6).')

                npts = int((nargs2-1) / 2)

                if(npts < 3):
                    raise Exception('Polygon specification requires an even number of values (>=6).')

                for j in range(1, nargs2):
                    if self._is_number(args2[j]) is False:
                        raise Exception('Polygon coordinate pairs must be numbers.')

                ra  = []
                dec = []

                for j in range(npts):
                    lon = float(args2[2 * j + 1])
                    lat = float(args2[2 * j + 2])

                    coord = SkyCoord(lon * u.degree, lat * u.degree, frame=coordsys)

                    raj  = coord.icrs.ra.deg
                    decj = coord.icrs.dec.deg

                    ra.append (raj)
                    dec.append(decj)


                retval = spt.polygon_search(npts, ra, dec, self.mode, self.level, self.xcol, self.ycol, self.zcol, self.indxcol, self.encoding)


                if(val == '0'):
                    retval['geom_constraint'] = retval['geom_constraint'].replace('<=', '>')
                    geomstr = "(" + retval['geom_constraint'] + ")"

                else:
                    geomstr = "(" + retval['geom_constraint'] + ") AND (" + retval['index_constraint'] + ")"


            # Point in box

            elif(func1 == 'point'
             and nargs1 == 3
             and self._is_number(args1[1]) is False
             and self._is_number(args1[2]) is False
             and func2 == 'box'):

                try:
                    coordsys = self._frame_lookup(args1[0])
                except Exception as e:
                    raise Exception(str(e) + ' in point() function.')

                if(nargs2 != 5):
                    raise Exception('box specification requires four coordinate pairs.')

                try:
                    coordsys = self._frame_lookup(args2[0])
                except Exception as e:
                    raise Exception(str(e) + ' in box() function.')

                for j in range(1, nargs2):
                    if self._is_number(args2[j]) is False:
                        raise Exception('Box parameters must be numbers.')

                lon = float(args2[1])
                lat = float(args2[2])

                width  = float(args2[3])
                height = float(args2[4])

                box = self._box(lon, lat, width, height)

                ra  = []
                dec = []

                for j in range(4):
                    lon = box[2 * j]
                    lat = box[2 * j + 1]

                    coord = SkyCoord(lon * u.degree, lat * u.degree, frame=coordsys)

                    raj  = coord.icrs.ra.deg
                    decj = coord.icrs.dec.deg

                    ra.append (raj)
                    dec.append(decj)

                npts = 4

                retval = spt.polygon_search(npts, ra, dec, self.mode, self.level, self.xcol, self.ycol, self.zcol, self.indxcol, self.encoding)

                if(val == '0'):
                    retval['geom_constraint'] = retval['geom_constraint'].replace('<=', '>')
                    geomstr = "(" + retval['geom_constraint'] + ")"

                else:
                    geomstr = "(" + retval['geom_constraint'] + ") AND (" + retval['index_constraint'] + ")"

            else:
                raise Exception('Invalid CONTAINS() clause')


        # DISTANCE FUNCTION

        elif funcName == 'distance':

            func1 = self.funcData[i]['args'][0]['name']
            func2 = self.funcData[i]['args'][1]['name']

            args1 = self.funcData[i]['args'][0]['args']
            args2 = self.funcData[i]['args'][1]['args']

            if func1 != 'point' or func1 != 'point':
                raise Exception('DISTANCE function arguments must be POINTs.')

            nargs1 = len(args1)
            nargs2 = len(args2)

            print(func1, args1, nargs1)
            print(func2, args2, nargs2)

            geomstr = ''


            # Database column first

            if(     func1 == 'point'
                and nargs1 == 3
                and self._is_number(args1[1]) is False
                and self._is_number(args1[2]) is False

                and func2 == 'point'
                and nargs2 == 3
                and self._is_number(args2[1]) is True
                and self._is_number(args2[2]) is True ):

                coordsys = ''
                try:
                    coordsys = self._frame_lookup(args1[0])
                except Exception as e:
                    raise Exception(str(e) + ' in point() function.\n')

                racol  = args1[1]
                deccol = args1[2]

                if(racol != 'ra'):
                    raise Exception('Invalid ra column reference in point() function.\n')

                if(deccol != 'dec'):
                    raise Exception('Invalid dec column reference in point() function.\n')

                try:
                    coordsys = self._frame_lookup(args2[0])
                except Exception as e:
                    raise Exception(str(e) + ' in point function.\n')

                lon = float(args2[1])
                lat = float(args2[2])

                if(lat < -90.):
                    raise Exception('Invalid latitude value in circle() function.')

                if(lat > 90.):
                    raise Exception('Invalid latitude value in circle() function.')

                try:
                    coord = SkyCoord(lon * u.degree, lat * u.degree, frame=coordsys)
                except Exception as e:
                    raise Exception(str(e))

                ra  = coord.icrs.ra.deg
                dec = coord.icrs.dec.deg

                try:
                    dtr = math.atan(1.) / 45.

                    x = math.cos(ra * dtr) * math.cos(dec * dtr)
                    y = math.sin(ra * dtr) * math.cos(dec * dtr)
                    z =                      math.sin(dec*dtr)

                    geomstr = 'acos(x*(' + str(x) + ')+y*(' + str(y) + ')+z*(' + str(z) + '))/' + str(dtr)

                    print(geomstr)

                except Exception as e:
                    raise Exception(str(e))


            # Database column second

            if( func1 == 'point'
            and nargs1 == 3
            and self._is_number(args1[1]) is True
            and self._is_number(args1[2]) is True

            and func2 == 'point'
            and nargs2 == 3
            and self._is_number(args2[1]) is False
            and self._is_number(args2[2]) is False ):

                coordsys = ''
                try:
                    coordsys = self._frame_lookup(args1[1])
                except Exception as e:
                    raise Exception(str(e) + ' in point() function.\n')

                racol  = args2[1]
                deccol = args2[2]

                if(racol != 'ra'):
                    raise Exception('Invalid ra column reference in point() function.\n')

                if(deccol != 'dec'):
                    raise Exception('Invalid dec column reference in point() function.\n')

                try:
                    coordsys = self._frame_lookup(args1[0])
                except Exception as e:
                    raise Exception(str(e) + ' in point function.\n')

                lon = float(args1[1])
                lat = float(args1[2])

                if(lat < -90.):
                    raise Exception('Invalid latitude value in circle() function.')

                if(lat > 90.):
                    raise Exception('Invalid latitude value in circle() function.')

                try:
                    coord = SkyCoord(lon * u.degree, lat * u.degree, frame=coordsys)
                except Exception as e:
                    raise Exception(str(e))

                ra  = coord.icrs.ra.deg
                dec = coord.icrs.dec.deg

                try:
                    dtr = math.atan(1.) / 45.

                    x = math.cos(ra * dtr) * math.cos(dec * dtr)
                    y = math.sin(ra * dtr) * math.cos(dec * dtr)
                    z =                    math.sin(dec * dtr)

                    geomstr = 'acos(x*(' + str(x) + ')+y*(' + str(y) + ')+z*(' + str(z) + '))/' + str(dtr)

                    print(geomstr)

                except Exception as e:
                    raise Exception(str(e))


        # Unrecognize function

        else:
            raise Exception('Invalid function: ' + funcName)

        return geomstr


    def _is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False


    def _frame_lookup(self, str_in):

        '''
        _frame_lookup() is a simple utility function to convert the
        coordinate system frame of reference as given in an ADQL
        geometric object specification to the equivalent string
        needed by the astropy SkyCoord() object initialization.

        Parameters
        ----------
        str_in : string, required
            ADQL coordsys string.

        Returns
        -------
        Coordsys string compatible with astropy SkyCoord().
        '''

        coordstr = str_in.strip("'")
        coordstr = coordstr.split(' ')[0].lower()

        retstr = ''
        for key in self.fconvert:

            if (coordstr.strip() == key.strip()):
                retstr = self.fconvert[key]
                break

        if (len(retstr) == 0):
            msg = 'Invalid coordinate system: ' + coordstr

            raise Exception(msg)

        return retstr


    def _box(self, lon, lat, width, height):

        '''
        _box() determines the polygon coordinates for the corners
        of an ADQL 'box'.

        Parameters
        ----------
        lon : double, required
        lat : double, required
            Coordinates of the box center,
        width : double, required
        height : double, required
            Great circle width and height of the box.

        Returns
        -------
        Array of four lon,lat pairs for the box corners.
        '''

        w = math.radians(width)
        h = math.radians(height)

        N = []
        V = []

        n = np.array([-math.sin(h / 2.), 0.,  math.cos(h / 2.)])
        N.append(n)

        n = np.array([-math.sin(w / 2.),  math.cos(w / 2.), 0.])
        N.append(n)

        n = np.array([-math.sin(h / 2.), 0., -math.cos(h / 2.)])
        N.append(n)

        n = np.array([-math.sin(w / 2.), -math.cos(w / 2.), 0.])
        N.append(n)

        for i in range(4):

            inext = (i + 1)%4

            v = np.cross(N[inext], N[i])

            v = v / np.linalg.norm(v)

            V.append(v)

        r = math.radians(lon)
        d = math.radians(lat)

        sr = math.sin(r)
        cr = math.cos(r)

        sd = math.sin(d)
        cd = math.cos(d)

        rot = np.array([[cr * cd,  -sr, -cr * sd],
                        [sr * cd,   cr, -sr * sd],
                        [     sd,   0.,       cd]])

        polygon = []

        for i in range(4):

            c = np.matmul(rot, V[i])

            lonc = math.degrees(math.atan2(c[1], c[0]))
            latc = math.degrees(math.asin(c[2]))

            polygon.append(lonc)
            polygon.append(latc)

        return polygon


    def sql(self, adql):

        """
        sql() converts an ADQL statement into SQL appropriate to an
        Oracle database with HTM or HPX spatial index columns.

        Parameters
        ----------
        adql : string, required
            ADQL string, including geometric constraints.

        Returns
        -------
        SQL string for Oracle, including pure-SQL spatial constraints.
        """


        # Look for a 'TOP <n>' construct in query.  Analyze it and
        # remove it from the query.  It will get added back in
        # later as part of (or all of) the query WHERE clause.

        patched_adql = adql

        patched_adql = re.sub(r'contains\s*\(', r'contains(', patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub(r'polygon\s*\(',  r'polygon(',  patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub(r'circle\s*\(',   r'circle(',   patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub(r'point\s*\(',    r'point(',    patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub(r'box\s*\(',      r'box(',      patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub(r'distance\s*\(', r'distance(', patched_adql, flags=re.IGNORECASE)

        tags = patched_adql.split(' ')

        ntags = len(tags)

        selectstr = ''
        topstr    = ''
        countstr  = ''
        indx_top = -1
        indx_count = -1

        for i in range (ntags):

            if (tags[i].lower() == 'select'):
                selectstr = tags[i].lower()

            if (tags[i].lower() == 'top'):
                topstr = tags[i].lower()
                indx_top = i

                l = i + 1
                while (len(tags[l]) == 0):
                    l = l + 1

                countstr = tags[l].lower()
                indx_count = l

        if(selectstr != 'select'):
            raise Exception('Query is not a SELECT statemement.')

        haveTop = False
        count   = 0

        if (topstr == 'top'):
            haveTop = True
            try:
                count = int(countstr)
            except ValueError:
                raise Exception('TOP counter is not an integer.')

        newadql = ''
        for i in range(ntags):

            if (i == indx_top or i == indx_count):
                continue

            newadql = newadql + tags[i]

            if(i < ntags-1):
                newadql = newadql + ' '


        # Parse (recursively) the input ADQL

        res = sqlparse.parse(newadql)

        stmt = res[0]

        self._parseADQL(stmt)


        # Check GEOM references for value.  Must look
        # both after and before ("GEOM = 1" or "1 = GEOM")
        # and blank out the right fields.

        ntoken = len(self.adql_tokens)
        ngeom = -1

        for i in range(len(self.adql_tokens)):

            if(self.adql_tokens[i] == 'GEOM'):

                ngeom = ngeom + 1


                # Check for "= 1" after 'GEOM'

                j = i + 1

                # Skip spaces

                while (j < ntoken):
                    if(self.adql_tokens[j] != ' '):
                        break

                    j = j + 1

                # Check for equals

                if(j < ntoken and self.adql_tokens[j] == '='):

                    j = j + 1

                    # Skip spaces again

                    while (j < ntoken):
                        if(self.adql_tokens[j] != ' '):
                            break

                        j = j + 1

                    # Check for 1/0

                    if(j < ntoken and (self.adql_tokens[j] == '0' or self.adql_tokens[j] == '1')):

                        # Remember this equality

                        self.funcData[ngeom]['val'] = self.adql_tokens[j]

                        # and blank out old representation

                        for k in range(i + 1, j + 1):
                            self.adql_tokens[k] = ''

                        continue


                # If we get this far, check for "1 =" before 'GEOM' instead

                j = i-1

                # Skip spaces

                while (j >= 0):
                    if(self.adql_tokens[j] != ' '):
                        break

                    j = j-1

                # Check for equals

                if(j >= 0 and self.adql_tokens[j] == '='):

                    j = j-1

                    # Skip spaces again

                    while (j >= 0):
                        if(self.adql_tokens[j] != ' '):
                            break

                        j = j-1

                    # Check for 1/0

                    if(j >= 0 and (self.adql_tokens[j] == '0'
                                   or self.adql_tokens[j] == '1')):

                        # Remember this equality

                        self.funcData[ngeom]['val'] = self.adql_tokens[j]

                        # and blank out old representation

                        for k in range(j, i):
                            self.adql_tokens[k] = ''

                        continue


        if(self.debug):
            self.debugfile.write('\n')
            self.debugfile.write('\n')
            self.debugfile.write('DECOMPOSED QUERY\n')
            self.debugfile.write('\n')
            self.debugfile.write('=========================================================================\n')
            self.debugfile.write('adql_tokens: \n')

            for i in range(len(self.adql_tokens)):
                self.debugfile.write('token ' + str(i) + ':   [' +  self.adql_tokens[i] + ']\n')

            self.debugfile.write('\n')
            self.debugfile.write('=========================================================================\n')

            self.debugfile.write('funcData:\n')
            self.debugfile.write(str(self.funcData))
            self.debugfile.write('\n')
            self.debugfile.write('=========================================================================\n')
            self.debugfile.write('\n')
            self.debugfile.write('\n')



        # Check for the location of the start and end of the WHERE clause.
        # Not needed if we don't have a 'TOP n' constraint.  An oddity:
        # The same sqlparse library under Linux and OSX parse differently;
        # e.g. in one "group by" is a single token and in the other it is three.

        indx = 0
        outstr = ''

        where_start = -1
        where_end   = -1

        imax = len(self.adql_tokens) - 1

        if(haveTop):

            for i in range(imax):

                if(where_start == -1 and self.adql_tokens[i].lower() == 'where'):
                    where_start = i

                if(where_end == -1 and i < imax
                and (   self.adql_tokens[i].lower() == 'group'
                     or self.adql_tokens[i].lower() == 'order')
                and self.adql_tokens[i + 1].lower() == ' '):

                    if(i + 2 < imax):
                        for j in range(i + 2,imax + 1):
                            if(self.adql_tokens[j].lower() == 'by'):
                                where_end = i
                                break
                            elif(self.adql_tokens[j].lower() != ' '):
                                break

                elif(where_end == -1 and i < imax
                and (   self.adql_tokens[i].lower() == 'group by'
                     or self.adql_tokens[i].lower() == 'order by')
                and self.adql_tokens[i + 1].lower() == ' '):
                    where_end = i

                elif(where_end == -1 and i < imax
                and  self.adql_tokens[i].lower() == 'having'):
                    where_end = i

                if(where_start == -1 and where_end == -1
                                     and self.adql_tokens[i + 1].lower() == ';'):
                    where_end = i


            if(where_start != -1 and where_end == -1):
                where_end = imax + 1

            if(self.debug):
                self.debugfile.write('\n')
                self.debugfile.write('WHERE clause: \n')
                self.debugfile.write('\n')
                self.debugfile.write('where_start: ' + str(where_start) + '\n')
                self.debugfile.write('where_end:   ' + str(where_end)   + '\n')
                self.debugfile.write('\n')
                self.debugfile.write('=========================================================================\n')


        # Create the output local SQL statement.  Where we have inserted a
        # placeholder GEOM tag, replace with a generated geometry/indexing
        # constraint.  If we have a 'TOP n' directive, add it to the WHERE
        # clause (or if there was no where clause, add one.

        for i in range(len(self.adql_tokens)):

            # GEOM placeholder processing

            if(self.adql_tokens[i] == 'GEOM'):

                constraint = self._geomConstraint(indx)

                outstr = outstr + constraint
                indx = indx + 1

            # Start of 'WHERE'

            elif(haveTop and i == where_start):
                outstr = outstr + self.adql_tokens[i] + ' ( '

            # End of 'WHERE'

            elif(haveTop and i == where_end):
                if(where_start == -1):
                    outstr = outstr + ' WHERE ROWNUM <= ' + str(count) + ' ' +  self.adql_tokens[i]
                else:
                    outstr = outstr + ' ) AND ROWNUM <= ' + str(count) + ' ' +  self.adql_tokens[i]

            # Everything else

            else:
                outstr = outstr + self.adql_tokens[i]

        if(haveTop and where_start == -1 and where_end == -1):
            outstr = outstr + ' WHERE ROWNUM <= ' + str(count)

        if(haveTop and where_start != -1 and where_end == imax+1):
            outstr = outstr + ') AND ROWNUM <= ' + str(count)

        if(self.debug):
            self.debugfile.write('\n')
            self.debugfile.write(f'output:  {outstr:s}\n')

        return outstr
