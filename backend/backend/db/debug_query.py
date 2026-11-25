from sqlalchemy.dialects import postgresql
from sqlalchemy.orm import Query

def stringify_query(query: Query) -> str:
    return str(query.statement.compile(
        dialect=postgresql.dialect(),
        compile_kwargs={"literal_binds": True}))