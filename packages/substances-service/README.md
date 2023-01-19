# Substances Microservice

## Stack

### Self-Managed

- FastAPI
- Postgresql
- gRPC
- GraphQL
- RabbitMQ
- 

### Serverless

- REST API for sync commands with using **AWS Lambda and API Gateway**
- GraphQL API for sync queries with using **AWS AppSync**
- Data persistence with using **AWS DynamoDB**
- Decouple microservices with events using **AWS EventBridge**
- Cloud stack development with IaC using **Serverless Framework**

## TODO

- Instructions for use AWS CloudFormation CDK instead of Serverless Framework.

## Project Structure

- `susbtances`
  - `core`
    - `substances_core (Poetry Package)`: Entities, Business Logic via Input/Output Ports defintion.
  - `database`: Data access adapters.
    - `dynamodb`
      - `substances_dynamodb (Poetry Package)`: repository and unit of work implementation for **AWS DynamoDB**.
    - `inmemory`
      - `substances_inmemory (Poetry Package)`: repository and unit of work implementation for a memory database (should be used for test only!)
    - `postgresql`
      - `susbtances_postgresql (Poetry Package)`:
  - `infrastructure`: Infrastructure adapters.
    - `etl`
      - `etl (Poetry Package)`:
    - `rabbitmq`
      - `substances_rabbitmq (Poetry Package)`:
    - `serverless`
      - `queue (Poetry Package)`: Message Queue implementation for **AWS SQS**.
  - `presentation`: Client interfaces.
    - `aws`
      - `api (Poetry Package)`: Serverless implementation for **AWS API Gateway**.
    - `graphql (Node.JS Package)`: Standalone Apollo GraphQL Server.
    - `grpc (Poetry Package)`: Standalone gRPC Server.
    - `rest (Poetry Package)`: FastAPI implementation.

## Developement

### Docker

**WARNING**: build context for any `docker-compose.yml` is `packages` path and should be taken into consideration when using `COPY` statement insides Dockerfile.

1. Start `docker-compose up`
2. Stop and remove `docker-compose down`

## Deployment

### AWS

#### Instructions

1. Move to AWS application folder `cd substances/presentation/aws`
2. Deploy with `serverless deploy`
