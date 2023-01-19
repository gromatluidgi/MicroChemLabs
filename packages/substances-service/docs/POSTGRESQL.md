# PostgreSQL

## Tools

- Docker
  - PostgreSQL
  - PGAdmin
- DBeaver

## References

- https://www.pgadmin.org/docs/pgadmin4/development/container_deployment.html

## Setup

### Docker

#### PostgreSQL

- ```docker pull postgres```
- ```docker run --name substances-postgresql -e POSTGRES_USER=admin -e POSTGRES_PASSWORD=admin -p 5432:5432 -v C:/_docker/_data/substances-postgresql:/var/lib/postgresql/data -d postgres```

Helper commands:
- docker start substances-postgresql
- docker stop substances-postgresql

#### PGAdmin

- ```docker pull dpage/pgadmin4:latest```

- ```docker run --name substances-pgadmin -p 5051:80 -e 'PGADMIN_DEFAULT_EMAIL=user@domain.com' -e 'PGADMIN_DEFAULT_PASSWORD=SuperSecret' -d dpage/pgadmin4```

- Navigate to http://localhost:5051 and logon with:
  - **Username**:user@domain.com
  - **Password**: SuperSecret

- Add a new server:
  - **Host**: substances-postgresql container address (```docker network inspect bridge```)
  - **Username**: admin
  - **Password**: admin

Helper commands:
- docker start substances-pgadmin
- docker stop substances-pgadmin


Notes:
- PGADMIN_DEFAULT_PASSWORD must at least be mixed with lower and uppercase; otherwise logon error could occur.