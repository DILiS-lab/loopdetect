from importlib.resources import files


def test_example_data_is_packaged():
    data_file = files("loopdetect.data").joinpath("li08_solution.tsv")
    assert data_file.is_file()