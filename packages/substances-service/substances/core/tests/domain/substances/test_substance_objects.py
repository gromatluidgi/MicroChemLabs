import pytest
from substances_core.domain.susbtances.objects import InchiKey


@pytest.mark.parametrize(
    "inchikey",
    [
        ("/ZOWIAZODTVDLOE"),
        ("QLNWXBAGRTU*KI-UHFFFAO_SA-@"),
        ("BQJCRHHNABKAKU-KBQPJGBKSA"),
        ("ZOWIAZODTVDLOE_GCRSARJPSA-N"),
    ],
)
def test_init_inchikey_raise_exception(inchikey):
    with pytest.raises(ValueError):
        InchiKey(inchikey)


@pytest.mark.parametrize(
    "inchikey",
    [
        ("QLNWXBAGRTUKKI-UHFFFAOYSA-N"),
        ("ZOWIAZODTVDLOE-GCRSARJPSA-N"),
    ],
)
def test_init_inchikey(inchikey):
    inchikey = InchiKey(inchikey)
    assert inchikey is not None
