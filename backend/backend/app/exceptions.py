"""
Backend Exceptions
"""


class GriException(Exception):
    pass


class SourceDBAlreadyExistsException(GriException):
    pass


class DomainAlreadyExistsException(GriException):
    pass


class SourceDBDoesNotExistException(GriException):
    pass


class MissingExplicitColourMapException(GriException):
    pass
