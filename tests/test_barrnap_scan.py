import pytest

from rdna_miner.steps.barrnap_scan import (
    barrnap_version_is_supported,
    require_barrnap_0_9,
)


def test_accepts_barrnap_0_9():
    assert barrnap_version_is_supported("barrnap 0.9\n")
    assert barrnap_version_is_supported("Barrnap 0.9")


def test_rejects_other_barrnap_releases():
    assert not barrnap_version_is_supported("barrnap 1.10.6")
    assert not barrnap_version_is_supported("barrnap 0.8")
    assert not barrnap_version_is_supported("barrnap 0.9.1")
    assert not barrnap_version_is_supported("0.9")
    assert not barrnap_version_is_supported("")


def test_require_barrnap_0_9_accepts_supported_version():
    require_barrnap_0_9("barrnap 0.9")


def test_require_barrnap_0_9_rejects_1_10():
    with pytest.raises(RuntimeError, match="0.9"):
        require_barrnap_0_9("barrnap 1.10.6")
