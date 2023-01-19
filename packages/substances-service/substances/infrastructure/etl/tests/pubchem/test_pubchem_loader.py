# pylint: disable=protected-access

import pytest
from shared.utils.path import PathUtils
from substances_core.domain.syncs.entities.item import SyncItem
from substances_etl.pubchem.loader import (
    PubchemSubstanceLoader,
    PubchemSubstanceLoaderOptions,
)


@pytest.fixture
def loader():
    data_dir = PathUtils.generate_temp_dir("MicroChemLabsTest")
    pubchem_loader = PubchemSubstanceLoader(
        PubchemSubstanceLoaderOptions(data_dir),
    )
    return pubchem_loader


def test_load(loader: PubchemSubstanceLoader):
    archive_path = PathUtils.join_path(
        PathUtils.generate_temp_dir("MicroChemLabsTest"),
        "Compound_050000001_050500000.sdf.gz",
    )
    item = SyncItem(
        token="65cd510f1d90c47b88e99b90b92de4e3",
        location=archive_path,
        name="Compound_050000001_050500000.sdf.gz",
        transformation=b"\x1f\x8b\x08\x00K\x94\xc1c\x02\xff\xed\x9d[s\x13\xc9\x15\xc7\xf3\x90\x0f2\xe5'\xbb\x96\x99L\xdf{R\xb5\x0ffX-\xb0F\x06\x93\x05\x02K\xb9\x84,@\x89,9\x92\xcd\xae\x93\xda\xaa|\x88\xe41\xdf5\xddsQ\x9f\x1e]\xb0\xfbt\xd9\xa6\x80b\x17\x1d\xcd\xcc\x99\xd3\xa7\xbb\xe7\xf2\xd3\x99\xff\xfc\xfb\x8f\xff=\xff\xdf\x1f\xaa?\xafN\xcd\x877\xff\xda9\xbb\xfc\xd3\xec\xdd\xdfF\xc3\xf3\x9d?';\x8b\x8bw\x8b\xf3\xc1t8Z\x1c\x0fg\xf3Qv2;\x1d\x8c\xa7\xd9\xe2b\xf1\xae\xfe:[\xae\x91=o?\xed\xdcKv\xea\x15\x8fG\x9fF\xd3\xf3\x85\xf1\xf4\xe6\xad\xf9r<\x1d~\x1c\x1f\xff}ti\xbe\xb8\xf6\x8e\xeau\x17\xd9#\xeb\xe4'\xe3\xc38\xac]\xed\xec\x1f\x96\xfdW\xcf_>\xd8\xbf\xff\xc3\xcb\xf20\xfd\xf9a\xaf\xd7\xdb?\xfc\xeb\xf3\xfd\xb4\xbf\xf3\xbbY\xed\xfdl~z1\x19`\xf6\xdak\\\x18o\xa7\xb3\xc9hh>\xcf\xad\x87\x92\x89\x87\xac('=\xda\x97\x87\xa2\xda\x9bY~\xfc\xebh\xfc\xe1c\xb5\x0bY\xa8\x8c\xda\xcd\x86\x83\xe9l:\x1e\x0e&\xc7\x8b\xd3\xf1d\xb4\xc0DS\xb6\xbe\x9e\xd7\xae\x8c\xfbO\x83\xc9\xc5\xa8\x8a\xa8\xec\xef\x96\xe5^\xb9\xfb\xfd\xe1\xde\x90\x0c\x87\xc3\xddr\xb27\x1c\x92~\xf5M\xd9\xfcC\xea\xd5\xec\xe7\xfe\x90\xda\x95\xea\x8fl\xb8\xdb3k\x0f\x87\xac\xfa\x87\x1aG%\xa9\x938\x9e~\x18\xcd\xcf\xe6c\xd0\xa1\x83\xc9\xf9h>\x1d\x9c\x8f?\x8d\x8e\xa7\x83S\xd7\xa8\xc5\xe8\xbcZ\xc7n8>1c`\xfc~<\x9a\xaf]\xbc\xb84m\xb9<]\xbbl8\x1f\r\xceG'\xc7\x83\xf3\xd5l\x9d\x98%\xe7\xe3S\x93\xa6\xe6\x83\xcd\xc2\xf1\xf1|tr1\x1c\x1d\x1f[\x1f\xd5\x16\xe7\x97g\xa3\xf5\xeb\x9b=\xbc\xd9\xb9\xff\xdd\xf0\xfe\x83\xa3\x1f\x06\x8f>\xa8\xef\xf6\x9f}\xff\xfd\xce\xdbj\xcf\x17g'n\xcf\xd3\x8b\xc9\xc46\xe4\xe2l0\xac\xda\x89\x1a\xbf\xd6K\x7fP\xc7\xbb\xec\xb5~\xfa\x86\xa6oD:\xfc8\x99\xcdg)MwO\xc6\xa3\xf3\x8f\x97\x93\xe1`\xfenp:\xbb\x9c\xec\r\xa6\xe3\xc9x:{k\x16\xce~\x9bU\x0b\xdf\xa6\xa4\xd9\xee\xcd.\xbd'\xd3\x93\xf1\xfb\xc9\x85\xd9\xfe\xec\xe3hj\xb6Xnk\xb79\xad\xb6X\xe7\xe4l|6\x9a\x8fO\xc6\xd3Q\xcaS\xbb\xcd\xec\xb7\xc1\xa9\xe9\xb5\xaa\xdb\xb1\xcd]\x1e\x15\xbe\x9c&\x9f\x0e\x16v4\xcaBfT2-(\xb7C\xfd\xd3\x87cof7\x83\xc2|GOlsD\x9es!\n\xf1\xcb4\xb1\x7f\x8e\x1e\xfc4>O\x96\x7f\xe8\x83_\xa6f\t/\x12A\x93$_\xfb\xb7(\x8a\xe4\x05\xcd\xf3\xbcv\xa1\xb3\\Rj?\xa9L\x19\xef\xf6S\x9e\x99\xc5yRN6\xf9\x80\x7f[7\x05\xd5\xd6M\xca\xbbnz\xc9\x95\xdd\x90<\x932g\xd6\x8dB\xb81M!\x85\xac\x1aE2\xea\xbb9\x0c\x8a\x86t\xa3\xb9\x86\x1b\x90b\x8eq\xc33Qh\x82vSd\xaa\xe0MO\xc5\xc8\x8d\xd9\xb8\xe3\xa6\x1f\x94\x1b\xdam\xd45\xdc\x10\x9a\xb1\x82\xae\xef\xa9\xebD#3\xc6\xf2\xf5)\xbe\x8e\x1b\x96)F\xf3\xdaaxn\xc0\xf0\x13\x98h\x96Ssu2\x94!\xe3&YIq\x19\x14M\x8e\x88\x06\x0c\xbf8\x8dJQ\xd1\x90L\x98c8\xd6\r\x98\x0c\x98\x14\x83hP\xb9q\x07Q\x86\xc9\x8d\x9b\x9a\x143n\\4\x1c\x13\x8dm\x8a\xae\xa2aq:<\x92\x1bL\xa3\xc0aK`\xc6\x8d\xcb\r*\xc5<#\x94sl\x87\x83q\xc31\xe3Fd\\r\xb2\xf6X\x1c6\xfcPn\xdc9\x1c\xd5S\xee\xcc\x10)7\n\xd3(\x9ai\x13\x0e\xb6Q\xe0\xac\xa9\xe2\xf4\x94\x8a3\xa7\x14\xa6Q\xb4\xde\x04;n\\\x8aQ\xd1\xb8\xd3\x9d\xc4\x9d\xa7\xda\x8b}\x11'\x1a\x15\xe5\xac\x19\xabQ\nw}\xd3F\x83h\x14Ixn\xff_[4\xe1\xc2Y,\xe1\xd2Y<\xa1\xc4\xaeQ[\"\xb1\xbbn-\x99P\xe1,\x950`\xe9\x84Ig\x15\t\xd1\xce\xa7\xb1\nhQ\x17\x0b\xc9\xab\xfd\x01\x8b9\x8bT{\x07\x16w\x16\xadb\x01\x96\xdb\x03aUd\xc0\xd2\x9e\x05\xd6\xe4U\xd4\xce\xe2`\x0f\"!\xd2\xb3\x14\xb4`\xd4\x12\xb6\xd6\xac\x06Zk\x9a\n\xda`\x1a\x07\xa26\xcd\xa1\xd2f|i)g\xc9\x84\xe5\x9e\xe5\xf6gV\xa3\xda-3\x99\xa7\x9e\xe5\xa26Ia\xcc-3\x16_Z\xc6=\xd8\x8ey9\xb3\x96\x8b\xc5l\xc4],f\x81=\xc4\xb7\x96N\xb8k\x9fI-w\xfdg6\x02k\x9a\xd4\xdaQ\x07,\xd7v\xb3\x80\xbb\xfd\x99\x05\xdc\xb5\xcf,\xe0\x85\xb34\xb0\x9e$\xc9\x0f\xfd\x07\xbfL\x1b\xbc\xc9N\x1c\xdfX\x8c\xfeq1\xb2P\xd7a\xb0\x0f\xd3\xd9|t<^,.*\xe0\xf7~0Y\x8c~\xbf\x97|\x91\xe8\xf8\xd9\xd1\xc3\xa3\x83\x97\xaf_\xdc\x7f\xd5/\x8f\x9e\xdf<:V\xe5\xa4'\xfa\xe2\x90\xafA\xc7\x8a\xd2\x8c\xdc(:nx\xf1\x92\nW0\xd8\xfeG\x01A\xde-\xf7\xfa\xa4,\xcbf\xbdr\xd7\xa3\xc8\xbd=\xf3\xb7a\xc7\x93j\xe5\xaf\x08!\xeb\x83\x8b\x0f\xb7\x89\x90+4\xda\xf0T\x91\xee\x9e\xcf\x1b8Z\x83P\x87SY\x83FSb\x91\xe8\xbb\x0b\xe3?\xa5i\x03WI\rW\xd9=~\x1d\xb8Zy:3\xab\xb6\xaen\x8f*\x7f1Yh@\xb3\xa2$\xa3\\p\xca\xe8\xd5@3\xc9\x85\xa6l+h\x16y\"X{\xe5B\xb6\x82fp\xaf+2FT(h^^\xde\x9b\xeb+\xed\xbb\xe9]\xeb2M4\x97i:#\x9a\x85\x83f\xb1\xbc\xf6\xe4\x82\x84G\xd3B1\x89\x89\xc6\xb9\xd1Q\xdc\xa4y7\xc5\x87Aw\x97\xb4\x1bM\x10=7\xf7\xba,<\x1a\x88vex4\xaeQy\xb7QA\xd8;\xa5\xddF\x85\x81f\x8a\x89\x06\xfe$\x80\x88\x06\xa2\x16D4n\x86\xaf\x0c\xbf\xa0\x1b\xf8\x94\xac\x1co\x82p\xc2J\x87\x07\xa1\x96H\x8d2\x1b\xcb\xf0h\x1c\x86B\xb9q\xa3\x18\x95b7\x8aQ)v\x87-\x8aI\xb1\xcb\r\xc1\xe4\xc6\xcd)\x8a\xc9\x8dk\x14\xc34\xcaE\x83\xe9)\xc0\x1dQ)v\x1d\xce0\x1d\xee\xa0\x18\x8b3\x19X\x9c\x0e\x17\x98\x0e\x87\xbf\xf2\xca\x18\xb9A\xbai\xf1%\xaa\xa7\\\x87\x8b83\\ \x86\x9f\xf7\xbb\xaa\x8e15%\xc6\x8d;kbr\x03\x1a%\xe3\x1c(\x14\xee\xb0\xe5((\x8bq\xba\x93\x98Q\x0c\x7f\xbcA\x0c?\xf8\x93@\x8c\xb3f\xacFi\xdc\x9cr\xf5?\xc1\xd1\x10Hi-\xaf\xe6\x1e\xaf\xe6\x90WCKT\xe4\x90,y\xb5\x00\xd4[y,[{,\xbb\xa8\xd8+]Rh\xee\x188!\x90\xe7Z\n\x9d{\x16\xf1(4\xf5,\x0e\xb93\xf5\x99t\x0e\xb93U\x1e\x93f\x904s\xe2Y\xc2\xe3\xcej\x13\x85\x96^,\x1e\x936\xab\x816\x98\x06\x00\xe6n\xad\x15BM<BM\x96L\x9aj\xcfr{\xb0l\x99l\"\xcd\xa4\"\xcd\xc0r{0\xabAB]\x93{`\x15\x90Iw\x08\xb5\xcb\x92q\x0f\x97i\x8f^\x17\x15\xaf\xe6\x8eP3\xcfr=fV\x83\xcbD0\xaf\x16n\xef\xb6\x8e.\xffF\xaf\x9f\xff\\\xbe~\xfd\xe2\xf5\xe3\xde\xd3\x83G\x07\x7f\xbdaz-\x1frZN6\xb0k\xc9\xf9\r\x97=[v\xdd\xa0hZ\xd1\xe9\xbd\xb2\xb4\xd8z\x13\xb7\xf6\xab\x9e\xab\x92\xe7%\xb5\xfejx\xf5\xab\xd1\xaf\xb7[\xf2\xec0mMW]\xedom\xb7,\xb6%\xad)Owy\xf3\x19\xb0UR\xb3\xd5)\xc4\xbbw\x88H\xdf\xa1v\xb6\xc5\xcd\x9ce\xb4\xa0\x82QqE\xe6\xac\x95$\x9f)n\xae/\x13\xae\xc2\x9c=\xb0\xd5\xde\xc8\x040g\x87\x92X\xd7M\x18;$\x99\x08w\xe3\xee\xc98\xc6\x8d#R\x02\xe3\x06\xf2\xb1\x8e\x9b~\xd0\xfd\xaa\xec\xa68\x8c\xab2\x8c\x1bw\xa5\x1f\xc7\x8d\xb9\x0bB\xb8q@te\xdc\x84\x02\xd1<\xcfc@\xbf<B4\xed\xc6\xd8\xdb\x17\x94\x1b\xf8\xc0\x00\"\xc5n\x14kL4\xb0\x12\x0e\x15M\x1b\x03\xca\r\xacf\x8f\x12\r\xca\x8d\x9b\x0c\x14\xe3\x06>w\x10e\x14\xa3\xdc\xb8qS\xc4q#\xe2\x0c?\x8eI\xb1;\xdd\xb1(=\xb5\xe6\x1c\x1e\x14\r\xaaQ\xb0\x807J\x8aY\x9chX\x9cQ,\xa2\x1cDS\x8e\xebp\xf7(Y\x94\x03\x85\x88\x13\x8d\x883\x8ae\x94\x9e\xc2\xb9q\xd7~*\xca\x0cOu\x94Q\xbcz\x99\x14\xe6\xa6\x8825\xd38gMLn<\xda\x96P\x8f\x88\xb2\x8aWRW\x07,a\x1d\xb0eh\xae\x0e\x98pHY\x89\xf0,\xe91W\x9f\xc0\x12\xcf\x02kj\xc8+\xad\xe5W\x0cs\xcf\xd2\xb0b\x98{\xd5\xc4\x90\x96\x12\x0ba\xa1\xb5\xb1\xb6\x98\xc2\x16\x99\xd5\x88GYaM\xb2\xf4\xa8\xae\xda\\\xf9k,H|=\xe6jV\x03\x99\xb0\xecT\xc1\xea^\xc0N\xad\xe5W\xf7\xe6\x9e\xe5W\xf7\n\xcf\x92\x1e\x81\xf5y,\xf7\x98\xab\xf3\xd2a\xa7\xd2#\xb0\xca#\xb0\nf\xde\xac\x06\x96\x99N\x80\x94\x95x$\xb5\xa6\xf9\xad\xc5\xb6W\x05\xcbo\xb4\xb4w\xffi\xef\xf0\xe0\xe1\x8f\x8f\x1f\x1f=9zp\xc3\xb4\xb4x\xc8e9\xe9mP\x89P\x8ce\xec\xb6pi\xb9,\xf5\x1dn\x03\xa6\xe5\xaa\\D\xc9\xcc\x97\xe5W\xc8M'\x8b[\xe5\xa6\xccG\x8a\xc3\xcb\xe1dvfR\xb7V;\x81\xb5|\xd0QF\x9e\xbe\xd9\xe5\xa9W\xd5z\xeaK/\xdc\xed\n\xdf/\xa0\xfdmm/\xa3\x99\xb9\xda\xa0\x94_QD\x82\x12N\xb8\xde^\xdbK\x13!\xaf\xcbY\x93\xd5\x8b\x9dI\xc8E\xf2\x8a\x9b\xb02X\x1a\x05\xd7\xae\x82\xae\xc3\x90\xa7j\x13\x0cg\x05\x0f\xbc\x0bL\xa3\x96\x88\x00G}!\xca\x0e'\x9b\xeeQ\xe1DE\xa1\xbe\xabw\xd3A\xd5\xb4)\xa6Q`2HL4\xeevF`\xdc8\xb2Ic\x01\xd2H\xbcL\xc4\x81\xc7Q\x08\x15&7f\x14\x17E\xce\xb1d\x93\x98y]\xd8R\x84\xfa@!D\xb0\x9b\"gE\x1dM\x91K\x15\xa1D\x18\xd1(\x99\xe5\x8a5h\x89\x16J\x06\xbb\x11\x8d\x1b\x99\x11\xc9d\x04\xea\x1b\x85\x88\x938\xc0\"\x8fB\xa8\x12\x0c \x05\xe7)\x19'\x1a\x0cKt\xea\x0f\x89\x8a\x03sPG?wH\xc7@@\xf7\x08\x14\xf6\x97\x94\xf6\x1c\x8e\x19\xc5 \xc5\x98\x1f\xab\xc0\x85I$zG\xe3\xa4\x18\x83$A\xa3t\x94\x1fd\xd28\x98\x1f\xe7f9\x05p\x93\xc1\xcd)\x1e\xc7\x8d\x883nd\x1c\x94\xad\xe2\x80u\x1d\x0b\xf3G\xf9\xc5_\xc7\xf9%E\xc5\xf9y\x08\x91\x9b\xae2\x86\xa0\xb0\xd2\x18\xd6\x08sX#l\x951\nH\xc4-Cu\xca\x18\xdcS\xc6\x00U\xba\xd6\xd2\x9b\xa8\xf7\x8aNF\xe1)c\x08\xcf\xd2\x1e\xcbV\xd0b\x1e\xd9\x86\xd5\xc4\xb4\xa2\xbb\x8ess\xaf&\x99{Z\x18\x90\x81sO\x0b\x83{\xb1\x88-U\xc8\xc5\x16\xea\x9dw\xeb\x8e\x1d\x936I\x81u\xc7\xd4\xe3\xe3\xcc[S\xc3JjK\xbd\x89W\x85\xec\xf6`+\x8dy\xb7\xee\x98;\x06.<Ko\xaa4\x16\x9e\xfa\x85\x80z\x17\x96\x96\xaf\xf0q_'\x03T!\xe7\x1e\x1f\xe7\x9b\xf88\x87}d+\x8d]&l\xa5\xf1\xd5\xeb\x8e]\xd4\xf6!k\x17\x8bY\xe0\xac\xaf\x97\xb2?}]>zy\xf8\xf4\xe0\xe9\xe3g?\x1e\xf6n\x98\xb2\xf3J\x8cy\x93\x143\xbfi)\xe6\xbeSa\x1e\xee\xb6\x02\xcct\rRg\xadNs\xb3\xd6\x90W\x15\xc9\xbc\xdc\x1b\xb2\xba\x8e\x99\x94_\tZ\xff\xf4\xfe\xd9-\xa3u\xde\xa2ef5\x86\x1b\x0c\xfcn4\xfdg\x85\x95O-\x0c\xde\x08\x97)(\xe1e\xe9\xee\xe9\x15\x04\x8do\x15\xa3\xdf\xa5\xb6.K\x93\xad\xee\xb2\xca\xb9$WC\xe6</rU|\xae4\xd9\x1e\xd3\xaf\xa4\xbb\x0c\xabx\x11\xba\xcb\xb0\xa83\x92(p\x14\x89b\x8c\xee\xb2'\xc5\x86p\x03\x9f\xac\x8c\x92\x1b\x8cD\xb1'k\x80\x10\x05\x86\x92\xa1\x88h \xfaD\xb8q\xb9\xd1q\x04\x93\x91J\xc7\xae49\x8a\x84\x1f\x8d\xa2t\x8c\x93\xe1\x85\n\x00Q\x94\x8eI\x1c\xe5Q\x16E\xa4\x13\x95\x1bO\x075\x8a\xd6g\x1e\xc7\r\x89#\n,\xe2\x88\xbbF\xeap\x19G\x05Z\xc6\xe9p\x15g2\x888\x12\xc5*N\x87\xa3\x94Ga10\x8d\xa3\xca\x1a%\x1a\xd4\xd4\x84\nGQ\xb4\x859N\xcd\xd7)\x8eDQ\xf3\xe5q\xce\x0cE\x1cM\xeaH\xd2\xd6\xa8\xa9\t5\xa9#\x89IG9\x88j\x8c\x1b\xa8\x1c\x11e2\xe88S\xb3\xc0H\x14\x034\x99\xd0\n\x85\x82\xd2\xe4\xdc\x03\xb1\x1c\x82X\xe6I\x14\xb3\x02\x82X\xbf\x18\x19\x82X\xd5\x01\xb1\x10\xbdj\xafP\xd9\x83\x98\xd6b\x1e\x88\xe5\xd0bd\x13\x88%\x9e\x0c1\xf5$\x8a\xa9'J\xcc\xb6\x14#3\xaf\xc4\x98o\x91(^\x11y\xa0\x1eP\xe5]Qb\xee, J,<Qb\xb1F\xf2\x01\x14*\xe7\xdd\xb2e\x80es\x0f\xbd2\xcf\xe2]\x19b\x80^\xe5\x96\xd2d\xd5\x05\xaa\xeb\x85\x87\xf3-\xc5\xc8\xa4+\xf9 \xbc\xd2\xe4o\xc5\xc8]\xe1\xe1\xbf\xf4\x9e<\xd9\x7fq\xf0\xecq\xf9\xe3\xc1O7_\x8c\xcc7\xbe\xb1N\x11\x9e\xe9\x9b\xaeE^GG\xeb\xe2b(\xd1\xd0bSn\x97\xf0%\x1b\xb55\xcb\xcdj\xa4)d\x1e\x92\xdeW\xc2K\xf9\xf9\x87\xdb\xe5\xa5\xeb\xcbo\xbb\xb2\x06o*\xc0xOD{a\x9b\xf5~v9\x9f\xcf&\xf5\xd7\xf4\xa6\xd9\xe9\x9dmw[zl\xe61c\x9cP\xae\xaf\xa8\xf1\xc0\xa9\xe4W\xa8=VW\xae=\xb6\xee\xea\x8b\"m\xae\xd6\x02\x8b\x86\xad2X\xae\xeb{L\x84\x1b\xab\xfdg\x7f\xf3\xb3`C\n\xafP\xedz\xe8\x92(\xca[\x1d\x03M\x83Eo\x9bF\xd9[\x0e\xbfQ\xd7\xc3\xba\x84\xd9+\x89Js\xb9(rt4\n\x13\r\xcb\x84\xa6\xebS|\xbd\xfaZ3lT\xcdG:\xd1\\\xcfMN\x18o\x10T!x0\x8f5mj\x1e\x17\xa4\xa2 \xc1n\x9a\xe1g\xe2b\xe1\x8dRYA\xe8\xfa\x9e\nJq\x9agD\x04Wm\xb2\x8c\xabV\xd8UrM\x82\xdd4\xe3\xc6l\xdc\xc9M\x19\x94b\x94\x9be\x95\x9b\x95R\x152\x1c%HY4\xbf\xffP.\x8b\xe0[\xc2\xb6\xa7\x08\xa6\xa7\\\x85\x19\xb15\xd19\xba\xa7H\x9c\x9eB\xb91\xf3\x9a\x92\x06y3sX\x0e\x8eF\x92\xa2y<\x99\xb2`7\xee\xb0\xb5rH\xbf\x9e\x1bU\x97\x98\xa7\xd5\xf93\x1clP\xdd\"\xa8B\xe7\xe1\xd1hb\xf9Hj\x01/\xd3\"8\xc5\xaa\x8a&]=\x87\x87\x8d\x1b\x86\x1b7JHQ\x9f\x19\n\xc9\xc3\x1b\x95k\xd6\x9c|E\x81 Y\xb4 \xba\x06\xf0\x92T\x9f\x02\xd9\x11\xd7\xb2&\xe7\x05\xa7*\xbc\xa7\xf2\xeaa\x00\x93\x1b)B\x9fu\xb0\xc3\x8f\xd3\xa6Q,\xa7:\xd8M{a\xc21\xe3\xc6\xbe\x02\"W\x8d\x024\xe5\x12\x1d\x8d\xc0D\xa3\xb2\x9cK\xbd\xf6R \xcc\x8d\x8cs,\x96\xb8h\xda\x0b\x13\x8e\xcbM\xeb&R\x8aU\x9cFiL4:S\xe6\xf0W\xdf\x12E\xe9\xf08n\x08\xea2i\xd9(\x9c\x9be\x8a1g\x06\xe2\xa9\x15SO\xad\xf83\xc8\x9bl\xa9=\x96^\xed\xb1\xd8\x02\xb9}\xacM7@nk\xc9-\xc8\xdb\x7fG\x9f\xf0\x90\xb7_{\xeci%\xfbo\xe5\x83\xb5\xc7\x0cV\xb9\xae\xd4\x1e\x13o\x19\xac/\x16^\xed\xb1\xf0\xe2\xec\xea!C\x88\xef\xe1\xfeN\x95\xb2\xf0j\x8f\xbb\x00\xbc\xf3V\xbebK]r\xee)n\xf8\xd5\xc6\xd2\x03\xe0@\x9b\x83wq\xb8\x0f\xc0\x8bM\xda\x1c\xba[{\xbc\x11\x87\xe7\xde\x1b\xfbV*\x91\x85g\xb9~\xe8\xbc\xcd\xef3\x95\xc8d\xbd\x02\xf2\xb7J\xe4u\x88\xfd\xc5\xb3\x83\xf2\xc5\xe3\xde\x8b\xa7G\xaf\x0f\x1e?\xbcy\xc4\x9eo\xd1\xfb\xa0\xea\x16\xe4\x91k4N\xaeY\x95\\\x7fn\n\x92\xf7\xea2\xe5\xaf\xab(\xf9\xd7\xfd\xfd\xef\xf6o\xb9(\x99\xb5\x85\xba|UBx[\x85\xae\xc5\xcf\xa2\x91\xbaX\x92f,\x83\xbe\xc1\x02\xe5\xbb\xda\xee\xe5\xbb\xfb\xcc\x9d\x8ebTqvE\xc6\xce\x94P\xf4\xea\xfa\x1e\x9f)V\xde\xf6t\\\x88\xbe\x07!\x18}\x0fX<\x88\x90\xc2\xd8\xa2\x1a\x11&\x13\xc2\xe3D\xa31\xd1\xb8\xc7!Qn\xb6\xe84\xf4\xe3H\xa8\x84iX\x08\x8c\x9b-\x0f\xae\x86\xe9{\xe8(\x8a\xd7)\xcaMt)\x0cr\x07dB\xe2\x08s\xc4\xd2>\x81j\xe0Q\x9e\xc6\xcd\xa3<\x8d\x9b\xc4\xd1\xc3\xc5=\x0c\x1e\xe7\xb1}p\xf4\x13w\xe0Ay\xf8\x86\xa6[\x97\xc2\x00n\xe2<~\x8ds\x13G\nc\x9bHRXn\xf4\x9d\xd2\"\x10q\xf4=\"\xa5\xb8\x883\xa7\"\x89\x96\xc88\xb2\xba*\x8ah\t\xc9\xa3LM\x94\x1b\x10\r\x893\x8aQ\x82.NB\xa5\x8825q\x8dZ\x8eb\x12\xe7\xcd\x0ci\x11CB%\xc5\x8d\x9b\xb6\xa7\x90n\\4$N4$\x86\xd8M\x8a\xe9\xa9nY9PSXe\xec\xc4\xd3\xf7\x10\x1ec\xa7\x1ecgW\xd5\xf7\xe8(zl\"\xee\x96\xaa\x93M\x8c\xfd3T\xddW\xfb\xf09:\xdbB\xd5\xb7\xe9{\xe4W\xd6\xf7 ^Y9\xeb\x96\x95\xfbo\x12\xe4\xde\x9b\x04\xb9\xf7&A\xbe\xeeM\x82\x9d\xb2r\xedick\xd8\xf6\x8e\xc65\xdb^V\xee)z@\xb5\x0f\xe9\x91s\xe9\xa9}t\x19;\xb4\n\x8f\xa3\xaf0v\xee1\xf6\x15\xaa\xbe\xfe-\x83wR\xdf\xe3\xed\x7f\xb2\xff\x03\x12\xde\xe9\x13\x7f\x9d\x00\x00",
    )
    result = loader.load(item)
    assert result is not None
    assert len(result) == 7


def test_load_from_pickle_file(loader: PubchemSubstanceLoader):
    archive_path = PathUtils.join_path(
        PathUtils.generate_temp_dir("MicroChemLabsTest"),
        "Compound_127500001_128000000.sdf.gz",
    )
    item = SyncItem(
        token="363c446d401b0ce2fa0a0316faac39a5",
        location=archive_path,
        name="Compound_127500001_128000000.sdf.gz",
        transformation="363c446d401b0ce2fa0a0316faac39a5.pickle",
    )
    result = loader.load(item)
    assert result is not None
