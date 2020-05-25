import sys
import math
import numpy as np
import logging
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

    
    def __init__(self, mode=SpatialIndex.HTM, level=7, xcol='x', ycol='y', \
        zcol='z', colname=None, encoding=None, debugfile=None):

        """
        ADQL() initialization allows you to set a number of parameters
        associated with the indexing, all of which have defaults.  

        Parameters
        ----------
        mode : integer, optional, default=SpatialIndex.HTM
            The spatial indexing supports both Heirarchical Triangular Mesh (HTM)
            and Hierarchical Equal Area isoLatitude Pixelization (HEALPix)
            tesselation of the sphere.
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
        colname : string, optional, default depends on mode and encoding
            though we recommend using a name based on mode and level, such
            as 'htm20' or 'hpx14'.
            Column name for spatial index.
        encoding : integer, optional, default depends on level
            Since the tesselations are a nesting of four values, one
            encoding is essentially "base 4" (e.g. 13112032) though the number
            is still stored as a decimal integer.  Early on (twenty years ago
            we stored the indices this way as it made debugging easier and we
            still support this encoding.  We would not recommend it for any
            new database. Use SpatialIndex.BASE10 instead.
        debugfile : string, optional, default None
            If set, defines the file where debug logging messages will go.
        """


        self.debug = False
        
        self.adql_tokens = []

        self.geomFuncs = []
        self.geomArgs  = []
        self.funcArgs  = []
        self.geomVals  = []

        self.inContains = 0
        self.inPoint    = 0
        self.inCircle   = 0
        self.inPolygon  = 0
        self.inBox      = 0

        self.mode     = mode
        self.level    = level
        self.xcol     = xcol
        self.ycol     = ycol
        self.zcol     = zcol


        self.colname  = colname
        self.encoding = encoding

        if(colname == None):

            if(encoding == None):
                self.encoding = SpatialIndex.BASE4
                self.colname = 'spt_ind'

            else:
                if(mode == SpatialIndex.HTM):
                    self.colname = 'htm' + str(level)
                else:
                    self.colname = 'hpx' + str(level)

        if(self.encoding == None):
            self.encoding = SpatialIndex.BASE10

        if (debugfile != None and debugfile != ''):
            self.debug = True
            self.debugfile = debugfile
            logging.basicConfig(filename=self.debugfile, filemode='w', level=logging.DEBUG)

        self.fconvert = {'icrs'          : 'icrs',
                         'fk5'           : 'fk5',
                         'fk4'           : 'fk4',
                         'j2000'         : 'fk5',
                         'b1950'         : 'fk4',
                         'ecliptic'      : 'geocentrictrueecliptic',
                         'galactic'      : 'galactic',
                         'supergalactic' : 'supergalactic'}


    def _parseADQL(self, stmt):

        # Pieces of ADQL statement.
        #
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
        # We also need to watch for 'TOP <n>' constructs.
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


        if self.debug:
            logging.debug('Enter __parseADQL')
            logging.debug('stmt:')
            logging.debug(stmt)
    
        for i in range(len(stmt.tokens)):
    
            token = stmt.tokens[i]
    
            if(self.debug):
                logging.debug('Token: [' + str(token) + ']')
        
    
            # STATEMENT
    
            # By the time we get here we should already be inside a
            # Statement but we'll keep this here for debugging purposes.
    
            if(isinstance(token, sqlparse.sql.Statement)):
                if(self.debug):
                    logging.debug('       TYPE: Statement')
                self._parseADQL(token)
        
    
    
            # COMMENT   
    
            # Our use case doesn't allow for comments but we'll leave this
            # in for debugging.
    
            elif(isinstance(token, sqlparse.sql.Comment)):
                if(self.debug):
                    logging.debug('       TYPE: Comment')
                    logging.debug("       VALUE: '" + token.value + "'")
        
    
    
            # IDENTIFIER 
    
            # Spotting one of the "contains()" Identifiers starts
            # special processing.
    
            elif(isinstance(token, sqlparse.sql.Identifier)):
                if(self.debug):
                    logging.debug('       TYPE: Identifier')
                    logging.debug("       VALUE: '" + token.value + "'")
    
                if(token.value == 'contains'):
                    self.inContains = 1
                    self.inPoint    = 0
                    self.inCircle   = 0
                    self.inPolygon  = 0
                    self.inBox      = 0
    
                    if(self.debug):
                        logging.debug('       inContains -> 1')
    
                    self.adql_tokens.append('GEOM')
    
                    if(self.debug):
                        logging.debug('       [GEOM] added to adql_tokens[](1)')
    
                elif(self.inContains and token.value == 'point'):
                    self.geomFuncs.append(token.value)
                    self.funcArgs = []
    
                    if(self.debug):
                        logging.debug(\
			    '[' + token.value + '] added to geomFuncs[](1)')
                        logging.debug('funcArgs[] initialized (point)')
    
                    self.inPoint   = 1
                    self.inCircle  = 0
                    self.inPolygon = 0
                    self.inBox     = 0
    
                    if(self.debug):
                        logging.debug('       inPoint -> 1')
    
                elif(self.inContains and token.value == 'circle'):
                    self.geomFuncs.append(token.value)
                    self.funcArgs = []
    
                    if(self.debug):
                        logging.debug(\
			    '[' + token.value + '] added to geomFuncs[](2)')
                        logging.debug('       funcArgs[] initialized (circle)')
    
                    self.inPoint   = 0
                    self.inCircle  = 1
                    self.inPolygon = 0
                    self.inBox     = 0
    
                    if(self.debug):
                        logging.debug('       inCircle -> 1')
    
                elif(self.inContains and token.value == 'polygon'):
                    self.geomFuncs.append(token.value)
                    self.funcArgs = []
    
                    if(self.debug):
                        logging.debug(\
			    '[' + token.value + '] added to geomFuncs[](3)')
                        logging.debug('       funcArgs[] initialized (polygon)')
    
                    self.inPoint   = 0
                    self.inCircle  = 0
                    self.inPolygon = 1
                    self.inBox     = 0
    
                    if(self.debug):
                        logging.debug('       inPolygon -> 1')
    
                elif(self.inContains and token.value == 'box'):
                    self.geomFuncs.append(token.value)
                    self.funcArgs = []
    
                    if(self.debug):
                        logging.debug(\
			    '[' + token.value + '] added to geomFuncs[](4)')
                        logging.debug('       funcArgs[] initialized (box)')
    
                    self.inPoint   = 0
                    self.inCircle  = 0
                    self.inPolygon = 0
                    self.inBox     = 1
    
                    if(self.debug):
                        logging.debug('       inBox -> 1')
    
    
                elif(self.inPoint or self.inCircle or self.inPolygon or self.inBox):
                    if(token.value != '(' and token.value != ' ' and token.value != ','):
                       self.funcArgs.append(token.value)
    
                       if(self.debug):
                           logging.debug('       [' + token.value + '] added to funcArgs[](1)')
    
                else:
                    self.adql_tokens.append(token.value)
    
                    if(self.debug):
                        logging.debug(\
			    '[' + token.value + '] added to adql_tokens[](1)')
        
    
    
            # IDENTIFIERLIST
    
            # If we are in the middle of processing a "contains()" function,
            # the IdentifierList is probably the function argument list.
    
            elif(isinstance(token, sqlparse.sql.IdentifierList)):
                if(self.debug):
                    logging.debug('       TYPE: IdentifierList')
        
                list = token.tokens
        
                count = len(list)
        
                if(self.debug):
                    logging.debug('list:')
                    logging.debug(list)
                    logging.debug(f'count= {count:d}')
        
                for j in range(count):
                
                    if(self.debug):
                        logging.debug('list[j]:')
                        logging.debug(list[j])
        
                    if(isinstance(list[j], sqlparse.sql.TokenList)):
                        self._parseADQL(list[j])
    
                    elif(self.inPoint):
                        if(list[j].value != ' ' and list[j].value != ','):
                           self.funcArgs.append(list[j].value)
    
                           if(self.debug):
                               logging.debug('       [' + list[j].value + '] added to funcArgs[](2)')
    
                    elif(self.inCircle):
                        if(list[j].value != ' ' and list[j].value != ','):
                           self.funcArgs.append(list[j].value)
    
                           if(self.debug):
                               logging.debug('       [' + list[j].value + '] added to funcArgs[](3)')
    
                    elif(self.inPolygon):
                        if(list[j].value != ' ' and list[j].value != ','):
                           self.funcArgs.append(list[j].value)
    
                           if(self.debug):
                               logging.debug('       [' + list[j].value + '] added to funcArgs[](4)')
    
                    elif(self.inBox):
                        if(list[j].value != ' ' and list[j].value != ','):
                           self.funcArgs.append(list[j].value)
    
                           if(self.debug):
                               logging.debug('       [' + list[j].value + '] added to funcArgs[](5)')
    
                    elif (self.inContains == 0):
                        self.adql_tokens.append(list[j].value)
    
                        if(self.debug):
                            logging.debug('       [' + list[j].value + '] added to adql_tokens[](2)')
    
        
    
            # WHERE
    
            # The where clause is often a combination of "contains()" 
            # and standard SQL constrains but the "where" itself 
            # should just be passed along.
    
            elif(isinstance(token, sqlparse.sql.Where)):
                if(self.debug):
                    logging.debug('       TYPE: Where')
                    logging.debug(token)
                self._parseADQL(token)
        
    
    
            # PARENTHESIS
    
            # We won't want the parentheses inside the "contains()"
            # processing to be reproduced in the output query but like 
            # a lot of these we need to keep the tokens otherwise.     
            # We'll handle this by disallowing the saving of tokens
            # in the generic TOKEN processing below under the right
            # conditions.
    
            elif(isinstance(token, sqlparse.sql.Parenthesis)):
                if(self.debug):
                    logging.debug('       TYPE: Parenthesis')
                self._parseADQL(token)
        
    
    
            # COMPARISON
    
            # There shouldn't be any comparisons inside the "contains()"
            # sections (though one right after them) so we can just 
            # process these normally.
    
            elif(isinstance(token, sqlparse.sql.Comparison)):
                if(self.debug):
                    logging.debug('       TYPE: Comparison')
                self._parseADQL(token)
        
    
    
            # TOKENLIST
    
            # Like PARENTHESIS, we can let the TOKEN logic deal with
            # things that shouldn't be output because of "contains()".
    
            elif(isinstance(token, sqlparse.sql.TokenList)):
                if(self.debug):
                    logging.debug('       TYPE: TokenList (generic)')
                self._parseADQL(token)
        
    
    
            #TOKEN
    
            #If we are handling "contains()" and encounter a token, 
            # it should either be a closing parenthesis (ending some
            # part of our functions) or a token we don't need (like
            # whitespace or commas).
    
            elif(isinstance(token, sqlparse.sql.Token)):
                if(self.debug):
                    logging.debug('       TYPE: Token (generic)')
        
                try:
                    if(token.value == ')'):
                        if(self.inContains):
    
                            if(self.inPoint):
                                self.inPoint = 0
                                self.geomArgs.append(self.funcArgs)
    
                                if(self.debug):
                                    logging.debug('       [' + self.funcArgs + '] added to geomArgs[](1)')
    
                            elif(self.inCircle):
                                self.inCircle = 0
                                self.geomArgs.append(self.funcArgs)
    
                                if(self.debug):
                                    logging.debug('       [' + self.funcArgs + '] added to geomArgs[](2)')
    
                            elif(self.inPolygon):
                                self.inPolygon = 0
                                self.geomArgs.append(self.funcArgs)
    
                                if(self.debug):
                                    logging.debug('       [' + self.funcArgs + '] added to geomArgs[](3)')
    
                            elif(self.inBox):
                                self.inBox = 0
                                self.geomArgs.append(self.funcArgs)
    
                                if(self.debug):
                                    logging.debug('       [' + self.funcArgs + '] added to geomArgs[](4)')
    
                            else:
                                self.inContains = 0
                        else:
                            self.adql_tokens.append(token.value)
    
                            if(self.debug):
                                logging.debug('       [' + token.value + '] added to adql_tokens[](3)')
    
                    elif(self.inPoint or self.inCircle or self.inPolygon or self.inBox):
                        if(token.value != '(' and token.value != ' ' and token.value != ','):
                           self.funcArgs.append(token.value)
    
                           if(self.debug):
                               logging.debug('       [' + token.value + '] added to funcArgs[](5)')
    
                    elif (self.inContains == 0):
                        self.adql_tokens.append(token.value)
    
                        if(self.debug):
                            logging.debug('       [' + token.value + '] added to adql_tokens[](4)')
    
    
                except:
                    pass
        
    
    
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
    
        if self.debug:
            logging.debug(f'Enter __geomConstraint: i= {i:d}')
    
        func1 = self.geomFuncs[2*i]
        func2 = self.geomFuncs[2*i+1]
    
        args1 = self.geomArgs[2*i]
        args2 = self.geomArgs[2*i+1]
    
        val   = self.geomVals[i]
    
    
        # Point in circle
    
        spt = SpatialIndex()

        nargs1 = len(args1)
        nargs2 = len(args2)

        if self.debug:
            logging.debug('')
            logging.debug(f'func1= {func1:s}')
            logging.debug(f'func2= {func2:s}')
            logging.debug(f'nargs1= {nargs1:d}')
            logging.debug(f'nargs2= {nargs2:d}')
   
        geomstr = ''

        if( func1 == 'point'  
        and nargs1 == 3 
        and self._is_number(args1[1]) == False 
        and self._is_number(args1[2]) == False

        and func2 == 'circle' 
        and nargs2 == 4 
        and self._is_number(args2[1]) == True  
        and self._is_number(args2[2]) == True  
        and self._is_number(args2[3]) == True ):
            
            coordsys = ''
            try:           
                coordsys = self._frame_lookup(args1[0])
            except Exception as e:
               raise Exception (str(e) + ' in point() function.')
    
            if self.debug:
                logging.debug('')
                logging.debug(f'point:coordsys= : {coordsys:s}')
            

            racol  = args1[1]
            deccol = args1[2]
    
            if(racol != 'ra'):
               raise Exception('Invalid ra column reference in point() function.')
    
            if(deccol != 'dec'):
               raise Exception('Invalid dec column reference in point() function.')
    

            try:           
                coordsys = self._frame_lookup(args2[0])
            except Exception as e:
                raise Exception (str(e) + ' in circle function.')
    
            if self.debug:
                logging.debug('')
                logging.debug(f'circle:coordsys= : {coordsys:s}')
     

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
                coord = SkyCoord(lon*u.degree, lat*u.degree, frame=coordsys)
            except Exception as e:
            
                if self.debug:
                    logging.debug('')
                    logging.debug(f'SkyCoord exception: e= {str(e):s}')
     
                raise Exception (str(e))
    
            ra  = coord.icrs.ra.deg
            dec = coord.icrs.dec.deg

            try:
                retval = spt.cone_search(ra, dec, radius, self.mode, \
                    self.level, self.xcol, self.ycol, self.zcol, \
                    self.colname, self.encoding)

            except Exception as e:
            
                if self.debug:
                    logging.debug('')
                    logging.debug(f'spt.cone_search exception: e= {str(e):s}')
     
                raise Exception (str(e))
    
            if(self.geomVals[i] == '0'):
                retval['geom_constraint'] = retval['geom_constraint'].replace('<=', '>')
                geomstr = "(" + retval['geom_constraint'] + ")"

            else:
                geomstr = "(" + retval['geom_constraint'] + ") AND (" + retval['index_constraint'] + ")" 

        # Point in polygon

        elif(func1 == 'point' 
         and nargs1 == 3 
         and self._is_number(args1[1]) == False 
         and self._is_number(args1[2]) == False

#         and func2 == 'polygon' and nargs2 >= 7):
         and func2 == 'polygon'):
	 
            if (nargs2 < 7):
                raise Exception ('Polygon specification requires at least 6 numbers.')

            try:           
                coordsys = self._frame_lookup(args1[0])
            except Exception as e:
               raise Exception (str(e) + ' in point() function.')

            if self.debug:
                logging.debug('')
                logging.debug(f'point:coordsys= : {coordsys:s}')
     
            try:           
                coordsys = self._frame_lookup(args2[0])
            except Exception as e:
               raise Exception (str(e) + ' in polygon() function.')
    
            if self.debug:
                logging.debug('')
                logging.debug(f'polygon: coordsys= : {coordsys:s}')
     

            if(nargs2%2 == 0):
                raise Exception('Polygon specification requires an even number of values (>=6).')

            npts = int((nargs2-1) / 2)

            if(npts < 3):
                raise Exception('Polygon specification requires an even number of values (>=6).')

            for j in range(1,nargs2):
                if self._is_number(args2[j]) == False:
                    raise Exception('Polygon coordinate pairs must be numbers.')

            ra  = []
            dec = []

            for j in range(npts):
                lon = float(args2[2*j+1])
                lat = float(args2[2*j+2])

                coord = SkyCoord(lon*u.degree, lat*u.degree, frame=coordsys)

                raj  = coord.icrs.ra.deg
                decj = coord.icrs.dec.deg
    
                ra.append (raj)
                dec.append(decj)


            retval = spt.polygon_search(npts, ra, dec, self.mode, self.level, self.xcol, self.ycol, self.zcol, self.colname, self.encoding)


            if(self.geomVals[i] == '0'):
                retval['geom_constraint'] = retval['geom_constraint'].replace('<=', '>')
                geomstr = "(" + retval['geom_constraint'] + ")"

            else:
                geomstr = "(" + retval['geom_constraint'] + ") AND (" + retval['index_constraint'] + ")" 


        # Point in box

        elif(func1 == 'point' 
         and nargs1 == 3 
         and self._is_number(args1[1]) == False 
         and self._is_number(args1[2]) == False
         and func2 == 'box'):
            
            try:           
                coordsys = self._frame_lookup(args1[0])
            except Exception as e:
               raise Exception (str(e) + ' in point() function.')

            if self.debug:
                logging.debug('')
                logging.debug(f'point: coordsys= : {coordsys:s}')
     

            if(nargs2 != 5):
                raise Exception('box specification requires four coordinate pairs.')

            try:           
                coordsys = self._frame_lookup(args2[0])
            except Exception as e:
               raise Exception (str(e) + ' in box() function.')
    
            if self.debug:
                logging.debug('')
                logging.debug(f'box: coordsys= : {coordsys:s}')
     

            for j in range(1,nargs2):
                if self._is_number(args2[j]) == False:
                    raise Exception('Box parameters must be numbers.')

            lon = float(args2[1])
            lat = float(args2[2])

            width  = float(args2[3])
            height = float(args2[4])
            
            box = self._box(lon, lat, width, height)

            ra  = []
            dec = []

            for j in range(4):
                lon = box[2*j]
                lat = box[2*j+1]

                coord = SkyCoord(lon*u.degree, lat*u.degree, frame=coordsys)

                raj  = coord.icrs.ra.deg
                decj = coord.icrs.dec.deg
   
                if self.debug:
                    logging.debug(f'j= {j:d} raj= {raj:f} decj= {decj:f}')


                ra.append (raj)
                dec.append(decj)

            npts = 4

            retval = spt.polygon_search(npts, ra, dec, self.mode, self.level, self.xcol, self.ycol, self.zcol, self.colname, self.encoding)

            if(self.geomVals[i] == '0'):
                retval['geom_constraint'] = retval['geom_constraint'].replace('<=', '>')
                geomstr = "(" + retval['geom_constraint'] + ")"

            else:
                geomstr = "(" + retval['geom_constraint'] + ") AND (" + retval['index_constraint'] + ")" 

        else:
            raise Exception('Invalid CONTAINS() clause')

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

        if self.debug:
            logging.debug ('')
            logging.debug (f'Enter __frame_lookup: str_in= {str_in:s}')

        coordstr = str_in.strip("'")
        if self.debug:
            logging.debug ('')
            logging.debug (f'coordstrstr= {coordstr:s}')
        
        coordstr = coordstr.split(' ')[0].lower()
	    
        if self.debug:
            logging.debug ('')
            logging.debug (f'coordstr= {coordstr:s}')

        retstr = ''
        for key in self.fconvert:

            if self.debug:
                logging.debug ('')
                logging.debug (f'key= {key:s}')
                logging.debug (f'value= {self.fconvert[key]:s}')
                logging.debug (f'key.strip= [{key.strip():s}]')
                logging.debug (f'coordstr.strip= [{coordstr.strip():s}]')
       
            if (coordstr.strip() == key.strip()):
                retstr = self.fconvert[key]
                break
	        
        if self.debug:
            logging.debug ('')
            logging.debug (f'retstr= {retstr:s}')

        
        if (len(retstr) == 0):
            msg = 'Invalid coordinate system: ' + coordstr 
        
            if self.debug:
                logging.debug ('')
                logging.debug (f'msg= {msg:s}')

            raise Exception (msg)

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

        n = np.array([-math.sin(h/2.), 0.,  math.cos(h/2.)])
        N.append(n)

        n = np.array([-math.sin(w/2.),  math.cos(w/2.), 0.])
        N.append(n)

        n = np.array([-math.sin(h/2.), 0., -math.cos(h/2.)])
        N.append(n)

        n = np.array([-math.sin(w/2.), -math.cos(w/2.), 0.])
        N.append(n)

        for i in range(4):

            inext = (i+1)%4

            v = np.cross(N[inext], N[i])

            v = v/np.linalg.norm(v)

            V.append(v)

        r = math.radians(lon)
        d = math.radians(lat)

        sr = math.sin(r)
        cr = math.cos(r)

        sd = math.sin(d)
        cd = math.cos(d)

        rot = np.array([[ cr*cd,  -sr, -cr*sd],
                        [ sr*cd,   cr, -sr*sd],
                        [    sd,   0.,     cd]])

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

        patched_adql = re.sub('contains\s*\(', 'contains(', patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub('polygon\s*\(',  'polygon(',  patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub('circle\s*\(',   'circle(',   patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub('point\s*\(',    'point(',    patched_adql, flags=re.IGNORECASE)
        patched_adql = re.sub('box\s*\(',      'box(',      patched_adql, flags=re.IGNORECASE)

        tags = patched_adql.split(' ')

        ntags = len(tags)

        if self.debug:
            logging.debug(f'ntags= {ntags:d}')
            
#            for i in range (ntags):
#                logging.debug (f'i= {i:d} tag={tags[i]:s}')


        selectstr = ''
        topstr    = ''
        countstr  = ''
        indx_top = -1
        indx_count = -1

        for i in range (ntags):
		
            if self.debug:
                logging.debug (f'i= {i:d} tag={tags[i]:s}')
        
            if (tags[i].lower() == 'select'):
                selectstr = tags[i].lower()

            if (tags[i].lower() == 'top'):
                topstr = tags[i].lower()
                indx_top = i

                l = i + 1
                while (len(tags[l]) == 0):
                    l = l + 1
          
                if self.debug:
                    logging.debug (f'l= {l:d} tag={tags[l]:s}')
        
                countstr = tags[l].lower()
                indx_count = l


        if self.debug:
            logging.debug (f'selectstr= {selectstr:s}')
            logging.debug (f'topstr= {topstr:s} indx_top={indx_top:d}')
            logging.debug (f'countstr= {countstr:s} indx_count={indx_count:d}')
        
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
        
        for i in range(len(self.adql_tokens)):
        
            if(self.adql_tokens[i] == 'GEOM'):
        
                # Check for "= 1" after 'GEOM'
        
                j = i+1
        
                # Skip spaces
        
                while (j < ntoken):
                    if(self.adql_tokens[j] != ' '):
                        break
        
                    j = j+1
        
                # Check for equals
        
                if(j < ntoken and self.adql_tokens[j] == '='):
        
                    j = j+1
        
                    # Skip spaces again
        
                    while (j < ntoken):
                        if(self.adql_tokens[j] != ' '):
                            break
        
                        j = j+1
        
                    # Check for 1/0
        
                    if(j < ntoken and (self.adql_tokens[j] == '0' or self.adql_tokens[j] == '1')):
        
                        # Remember this equality
        
                        self.geomVals.append(self.adql_tokens[j])
        
                        # and blank out old representation
        
                        for k in range(i+1, j+1):
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
        
                        self.geomVals.append(self.adql_tokens[j])
        
                        # and blank out old representation
        
                        for k in range(j, i):
                            self.adql_tokens[k] = ''
        
                        continue
        
        
        if(self.debug):
            logging.debug('============================')
            logging.debug('geomFuncs: ')
            logging.debug(self.geomFuncs)
            logging.debug('geomArgs: ')
            logging.debug(self.geomArgs)
            logging.debug('')
            logging.debug('============================')
            logging.debug('')
            logging.debug('geomVals: ')
            logging.debug(self.geomVals)
            logging.debug('')
            logging.debug('============================')
        outstr = ''
        
        indx = 0
        

        # Check for the location of the start and end of the WHERE clause.
        # Not needed if we don't have a 'TOP n' constraint.

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
                and self.adql_tokens[i+1].lower() == ' '):

                    if(i+2 < imax):
                        for j in range(i+2,imax+1):
                            if(self.adql_tokens[j].lower() == 'by'):
                                where_end = i
                                break;
                            elif(self.adql_tokens[j].lower() != ' '):
                                break;

                elif(where_end == -1 and i < imax 
                and  self.adql_tokens[i].lower() == 'having'):
                    where_end = i

                if(where_start == -1 and where_end == -1 
                                     and self.adql_tokens[i+1].lower() == ';'):
                    where_end = i


            if(where_start != -1 and where_end == -1):
                where_end = imax+1


        # Create the output local SQL statement.  Where we have inserted a
        # placeholder GEOM tag, replace with a generated geometry/indexing
        # constraint.  If we have a 'TOP n' directive, add it to the WHERE
        # clause (or if there was no where clause, add one.

        for i in range(len(self.adql_tokens)):

            if(self.debug):
                logging.debug('processing token: ' + self.adql_tokens[i])


            # GEOM placeholder processing

            if(self.adql_tokens[i] == 'GEOM'):

                if(self.debug):
                    logging.debug(f'call __geomConstraint: indx= {indx:d}')

                constraint = self._geomConstraint(indx)

                if(self.debug):
                    logging.debug('returned __geomConstraint')

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

            if(self.debug):
                logging.debug('')
                logging.debug(f'output:  {outstr:s}')


        if(haveTop and where_start == -1 and where_end == -1):
                outstr = outstr + ' WHERE ROWNUM <= ' + str(count)

        if(haveTop and where_start != -1 and where_end == imax+1):
                outstr = outstr + ') AND ROWNUM <= ' + str(count)

        
        if(self.debug):
            logging.debug('')
            logging.debug(f'output:  {outstr:s}')


        return outstr
