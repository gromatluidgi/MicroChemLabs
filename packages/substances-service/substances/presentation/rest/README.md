# Substances Service - Rest API

## Requirements

- Running instance of Postgresql Service
- Running instance of RabbitMQ Service

## Docker

- Build context: `../../../../`
- Container name: substances-rest-api
- Hostname: 
- Network: microchemlabs
- Environment variables:
  - DB_URI
  - RABBITMQ_HOST

### Instructions

1. Create network `docker network create microchemlabs` if not exists
2. Build `docker build -t substances-rest-api ../../../../ --file Dockerfile`
3. Run `docker run -d --name mcl-substances-rest-api --hostname mcl-substances-rest-api --network microchemlabs -p 5400:5400 --env RABBITMQ_HOST=172.19.0.2 --env DB_URI="postgresql://172.19.0.3/substances?user=admin&password=admin" substances-rest-api`
4. **Expected logs**:
```
2023-01-04 08:31:07 Skipping virtualenv creation, as specified in config file.
2023-01-04 08:31:08 INFO:     Started server process [1]
2023-01-04 08:31:08 INFO:     Waiting for application startup.
2023-01-04 08:31:08 INFO:     Application startup complete.
2023-01-04 08:31:08 INFO:     Uvicorn running on http://0.0.0.0:5400 (Press CTRL+C to quit)
2023-01-04 08:31:08 INFO:pika.adapters.utils.connection_workflow:Pika version 1.3.1 connecting to ('172.19.0.2', 5672)
2023-01-04 08:31:08 INFO:pika.adapters.utils.io_services_utils:Socket connected: <socket.socket fd=23, family=AddressFamily.AF_INET, type=SocketKind.SOCK_STREAM, proto=6, laddr=('172.19.0.4', 41004), raddr=('172.19.0.2', 5672)>
2023-01-04 08:31:08 INFO:pika.adapters.utils.connection_workflow:Streaming transport linked up: (<pika.adapters.utils.io_services_utils._AsyncPlaintextTransport object at 0x7fc7bc1344c0>, _StreamingProtocolShim: <SelectConnection PROTOCOL transport=<pika.adapters.utils.io_services_utils._AsyncPlaintextTransport object at 0x7fc7bc1344c0> params=<ConnectionParameters host=172.19.0.2 port=5672 virtual_host=/ ssl=False>>).
2023-01-04 08:31:08 INFO:root:Starting server on [::]:50051
2023-01-04 08:31:08 INFO:substances_rest.queues.consumer:Connection opened
2023-01-04 08:31:08 INFO:pika.adapters.utils.connection_workflow:AMQPConnector - reporting success: <SelectConnection OPEN transport=<pika.adapters.utils.io_services_utils._AsyncPlaintextTransport object at 0x7fc7bc1344c0> params=<ConnectionParameters host=172.19.0.2 port=5672 virtual_host=/ ssl=False>>
2023-01-04 08:31:08 INFO:pika.adapters.utils.connection_workflow:AMQPConnectionWorkflow - reporting success: <SelectConnection OPEN transport=<pika.adapters.utils.io_services_utils._AsyncPlaintextTransport object at 0x7fc7bc1344c0> params=<ConnectionParameters host=172.19.0.2 port=5672 virtual_host=/ ssl=False>>
2023-01-04 08:31:08 INFO:substances_rest.queues.consumer:Channel opened
```