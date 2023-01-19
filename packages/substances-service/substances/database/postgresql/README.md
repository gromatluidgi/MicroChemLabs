## Docker

### Dump using pg_dump

- Run `docker exec -i mcl-susbstances-postgresql /bin/bash -c "PGPASSWORD=admin pg_dump --username admin substances" > C:/_docker/_data/substances-postgresql/mcl_substances_dump.sql`

- Run `docker exec -i substances-postgresql /bin/bash -c "PGPASSWORD=admin pg_dump --username admin substances" > C:/_docker/_data/substances-postgresql/mcl_substances_dump.sql`

### Restore using psql

- Run docker exec -i pg_container_name /bin/bash -c "PGPASSWORD=pg_password psql --username pg_username database_name" < /path/on/your/machine/dump.sql
