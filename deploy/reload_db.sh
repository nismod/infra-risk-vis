#
# Load terracotta metadata database to remote - avoids re-running the ingest
#
set -ex

dbname=$1

# Ensure PG* environment variables are set for *remote* database
# PGHOST=
# PGDATABASE=
# PGUSER=
# PGPASSWORD=

# Create remote database
createdb $dbname

# Run pgloader to copy database
# May be best installed from source: https://github.com/dimitri/pgloader
# Could switch to plain dump/restore - this utility is/was helpful in
# translating from MySQL/SQLite metadata to PostgreSQL.
pgloader \
    "pgsql://global_dev:password@localhost:5432/$dbname" \
    "pgsql://$PGUSER:$PGPASSWORD@$PGHOST/$dbname?sslmode=allow"
