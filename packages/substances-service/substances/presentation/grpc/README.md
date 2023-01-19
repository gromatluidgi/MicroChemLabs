## Docker

- Build context: `../../../../`
- Container name: 
- Hostname: 
- Network: 
- Environment variables:
  - DB_URI
  - RABBITMQ_HOST

### Instructions

1. Create network `docker network create virtual-lab`
2. Build `docker build ../../../../ -t substances-grpc`
3. Run `docker run -d --rm --name vlab-substances-grpc --hostname vlab-substances-grpc --network virtual-lab -p 4200:4200  substances-grpc`