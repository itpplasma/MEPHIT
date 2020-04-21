from magdif_fffi import magdif_config, magdif

# TODO: doesn't work yet
def test_readmesh():
    magdif_config.read_config()
    print("read config")
    magdif.magdif_init()
    assert True
