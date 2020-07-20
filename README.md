# ADQL

## ADQL to Local SQL Translator

Develop branch (please do any new development from here)

ADQL is a variant of SQL that is understood by no DBMSs, the main difference 
being support for sky spatial constraints.  

One approach to supporting ADQL is to extend the DBMS with these ADQL-specific
constructs.  The alternate approach taken here is to translate the ADQL into
local SQL before sending it to the DBMS.

Since we also want to make use of spatial indexing to speed up spatial searches,
we include in this translation additional spatial index column constraints.  
This requires reference to special added columns. 

The hope is that this code can be leveraged by others who either wish to add
spatial indexing to their DBMS the same way we have or would like a starting
place for developing code around their own DBMS.

This package is complementary to but independent of our implemetation of the
IVOA Table Access Protocol (TAP); a web service for managing remote database
requests for astronomy in a general way.  TAP manages requests/results and is
generally fed ADQL.  This ADQL translator is used by our TAP implementation
to turn ADQL into SQL the local DBMS can handle directly.

Since this is pure Python, the PyPI package can be built with a simple
"<i>python setup.py bdist_wheel</i>".
