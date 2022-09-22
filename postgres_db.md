# Postgres DB

Using a Postgres DB we can:

- create a DB (i.e. pop_genomic)
- in the DB we can have various schemas representing the groups (e.g. soranzo, blagoje ...)
- in each schema we have a set of fixed tables (samples, analysis, projects, individuals, additional ids, etc). Permission are then assigned to users only for the relevant schema(s) so they can only interact with a certain part of the data
- we can create one admin schema containing views that concatenate the fixed tables across all schemas for center-wide view of data.

In this way we would be able to restrict access for group members so they can see only group related data, but admins can do queries across all groups data so we can have a view of what's going on in the whole center

## Implement Postgres on the HPC

I've been able to run postgres on the HPC using singularity and the [official postgres docker image](https://hub.docker.com/_/postgres). See data and script in /group/soranzo/edoardo.giacopuzzi/postgres/.
Essentially, we need a postegres service running and then we can connect to it using the same image as client.

### 1.Init the DB with postgres super-user

`PGDATA` must point to an empty folder, the superuser name is set by `-U`.

```bash
singularity exec \
    -B /group/soranzo/edoardo.giacopuzzi/postgres/data:/data/postgres \
    --env POSTGRES_PASSWORD="18-Set-1991" \
    --env PGDATA=/data/postgres \
    /project/alfredo/singularity/postgres_latest.sif initdb \
        --no-locale -U postgres
```

### 2.Set authentication method

First, login as postgres super-user and assign a password to it using `ALTER ROLE postgres ENCRYPTED PASSWORD 'mypassword'`.

Then set authentication method by changing the default `trust` in the `pg_hba.conf` file in the postgres data folder. You can use `md5` to request password typing ot `peer` to request that system user matches a corresponding DB user. See [postgres docs on atuthentication](https://www.postgresql.org/docs/9.1/auth-methods.html).

### 3.Run the server service

`PGDATA` here point to a previously created postgres data folder created using initdb

```bash
singularity exec \
    -B /group/soranzo/edoardo.giacopuzzi/postgres/data:/data/postgres \
    -B /group/soranzo/edoardo.giacopuzzi/postgres/var:/var/run/postgresql \
    --env PGDATA=/data/postgres \
    /project/alfredo/singularity/postgres_latest.sif \
    postgres -h 127.0.0.1 -p 9997
```

### 4.Connect to a running server

`-U` set the user name and `-d` the name of the database to connect to

```bash
user=$1
db_name=$2

singularity exec \
    /project/alfredo/singularity/postgres_latest.sif psql \
        -U ${user} -h 127.0.0.1 -p 9997 -d ${db_name}
```

## Manage permissions

We can combine users into groups and the assign specific permissions to the groups for the corresponding schema so all users in the group can read from the schema (representing research group data). We also have an admins group and schema for operations that need to traverse the whole DB.

### Alter password for a user

`ALTER USER xxx WITH ENCRYPTED PASSWORD 'yourpassword'`

### Create groups

`CREATE GROUP blagoje_group;`

### Create user in groups

`CREATE USER blagoje WITH IN GROUP blagoje_group ENCRYPTED PASSWORD 'blagoje';`

### Create schemas and tables in them

This creates a schema with a specific owner:
`CREATE SCHEMA admins AUTHORIZATION soranzo;`

Table in a schema:
`CREATE TABLE schema.table (val1 INTEGER, val2 TEXT)`

### Grant access to a schema

By default only the owner can access a new schema. To grant access to user or group

`GRANT USAGE ON SCHEMA myschema TO group/user;`

Note that this does not grant any access to tables present in the schema. Access to table data must also be set.

### Grant access to table

You can set permissions for a specific table for a user or group using

```bash
GRANT { { SELECT | INSERT | UPDATE | DELETE | TRUNCATE | REFERENCES | TRIGGER }
    [, ...] | ALL [ PRIVILEGES ] }
    ON { [ TABLE ] table_name [, ...]
         | ALL TABLES IN SCHEMA schema_name [, ...] }
    TO role_specification [, ...] [ WITH GRANT OPTION ]
    [ GRANTED BY role_specification ]
```

See the [GRANT docs](https://www.postgresql.org/docs/current/sql-grant.html) for full details.

## Tables

### Record time and user automatically

```postgres
CREATE TABLE mytable (
    create_time TIMESTAMP NOT NULL DEFAULT current_timestamp,
    create_by TEXT NOT NULL DEFAULT current_user
);
```

### Reference another table

```postgres
CREATE TABLE mytable (
    individual_id INTEGER NOT NULL REFERENCES individuals(id) ON UPDATE CASCADE
)
```

## Triggers

Triggers are action performed automatically after an INSERT, UPDATE or DELETE on a table. To use a trigger generally you need to first create a function that defines the operation to perform.

```postgres
CREATE OR REPLACE FUNCTION trigger_set_timestamp()
RETURNS TRIGGER AS $$
BEGIN
  NEW.updated_time = current_timestamp;
  NEW.updated_by = current_user;
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER set_updated_projects
BEFORE UPDATE ON projects
FOR EACH ROW
EXECUTE PROCEDURE trigger_set_timestamp();
```
