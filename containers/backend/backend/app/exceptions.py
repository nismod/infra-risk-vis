"""
Backend Exceptions
"""


class GriException(Exception):
    """Base"""

    pass


class SourceDBAlreadyExistsException(GriException):
    pass

class SourceDBDoesNotExistException(GriException):
    pass
