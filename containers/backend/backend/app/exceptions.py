"""
Backend Exceptions
"""


class GriException(Exception):
    """Base"""

    pass


class SourceDBAlreadyExistsException(GriException):
    pass


class DomainAlreadyExistsException(GriException):
    pass


class SourceDBDoesNotExistException(GriException):
    pass


class MissingExplicitColourMapException(GriException):
    pass
