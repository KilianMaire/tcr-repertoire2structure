from rep2struct import webapp


def test_parse_args_defaults():
    base, port, open_browser = webapp.parse_args([])
    assert base == "runs"
    assert port == 8000
    assert open_browser is True


def test_parse_args_base_and_port():
    base, port, open_browser = webapp.parse_args(["mydir", "--port", "9123"])
    assert base == "mydir"
    assert port == 9123
    assert open_browser is True


def test_parse_args_no_browser_flag():
    base, port, open_browser = webapp.parse_args(["--no-browser"])
    assert open_browser is False
    assert base == "runs"
