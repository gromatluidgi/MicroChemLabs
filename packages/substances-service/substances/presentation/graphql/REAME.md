# Substances - GraphQL Apollo Subgraph

- Federation quickstart (local subgraphs): https://www.apollographql.com/docs/federation/quickstart/local-subgraphs

## Docker

- Build context: `../../../../`
- Container name: 
- Hostname: 
- Network: 
- Environment variables:

1. Build `docker build ../../../../ -t substances-subgraph`
2. Run `docker run --rm -d --name vlab-substances-subgraph --hostname substances-subgraph --network virtual-lab -p 3100:4000 substances-subgraph`