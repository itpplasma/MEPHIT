from magdif_fffi import magdif_config, magdif

# TODO: doesn't work yet
def test_answer():
    magdif_config.log_open()
    magdif.read_mesh()
    assert True
