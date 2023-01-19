from rdkit.RDPaths import *
from _typeshed import Incomplete
from pyPgSQL import PgSQL as PgSQL
from pysqlite2 import dbapi2 as dbapi2

RDBaseDir: Incomplete
RDCodeDir: Incomplete
RDDataDir: Incomplete
RDDocsDir: Incomplete
RDDemoDir: Incomplete
RDBinDir: Incomplete
RDProjDir: Incomplete
RDContribDir: Incomplete
rpcTestPort: int
pythonTestCommand: str
defaultDBUser: str
defaultDBPassword: str

class ObsoleteCodeError(Exception): ...
class UnimplementedCodeError(Exception): ...

pythonExe: Incomplete
usePgSQL: bool
useSqlLite: bool
RDTestDatabase: str
RDDataDatabase: str
molViewer: Incomplete
