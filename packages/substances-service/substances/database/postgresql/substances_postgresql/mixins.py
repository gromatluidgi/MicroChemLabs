from sqlalchemy import Column, DateTime
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.orm import declarative_mixin
from sqlalchemy.sql import expression


class utcnow(expression.FunctionElement):
    """utcnow"""

    type = DateTime()
    inherit_cache = True


@compiles(utcnow, "postgresql")
def pg_utcnow(_, __, **___):
    return "TIMEZONE('utc', CURRENT_TIMESTAMP)"


@compiles(utcnow, "mssql")
def ms_utcnow(_, __, **___):
    return "GETUTCDATE()"


@declarative_mixin
class TimestampMixin(object):

    created_at = Column(DateTime, default=utcnow(), server_default=utcnow())
    updated_at = Column(
        DateTime, onupdate=utcnow(), server_default=utcnow(), server_onupdate=utcnow()
    )
