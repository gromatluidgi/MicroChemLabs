# Substances Service: DynamoDB

## Access Patterns

Access patterns or query patterns define how the users and the system access the data to satisfy business needs. 
[Get more infos](https://docs.aws.amazon.com/prescriptive-guidance/latest/dynamodb-data-modeling/step3.html).

| Access Patterns                                |
| :--------------------------------------------- |
| Get sync for a given syncId                    |
| Get substance for a given substanceId          |
| Get similar substances for a given substanceId |

## Tests

Moto: http://docs.getmoto.org/en/latest/docs/services/dynamodb.html?highlight=dynamodb

## Setup

- aws --profile=via dynamodb create-table --cli-input-json file://json/Syncs.json


## Coverage

poetry run pytest --cov