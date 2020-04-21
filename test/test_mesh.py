from magdif_fffi import magdif_config, magdif

# TODO: write tests
# TODO: read config and initialize with fixtures

def test_magdif_init():
    magdif_config.read_config()
    print("read config")
    magdif.magdif_init()
    assert True  # TODO: check something actually
