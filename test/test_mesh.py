import numpy as np
import pytest
from magdif_fffi import magdif_config, magdif

# TODO: write tests
# TODO: read config and initialize with fixtures


@pytest.fixture(scope='module', autouse=True)
def setup():
    magdif_config.read_config()
    print("Finished reading config")
    magdif.magdif_init()
    print("Magnetic field initialized")

def test_nan():
    assert not np.isnan(magdif.b0flux).any()
    assert not np.isnan(magdif.bnflux).any()
    assert not np.isnan(magdif.bnphi).any()
