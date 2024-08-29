def pytest_addoption(parser):
    parser.addoption(
        "--local",
        action="store_true",
        dest="local",
        default=False,
        help="Run local without docker.",
    )
