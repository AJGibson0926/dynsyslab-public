# DynSysLab Codex MCP harness

This is DynSysLab's local, personalized MCP development surface. It uses the
machine-wide official Python MCP SDK and defaults to stdio, so it is usable by
a local MCP client without exposing a listening network service.

## Safety boundary

The initial server supplies read-only orientation tools only. It never exposes
access keys, account metadata, arbitrary shell execution, archive mutation, or
external cloud controls. Future write-capable tools must be explicit,
capability-scoped, and delegate archive operations through Archiver.

## Run locally

```powershell
python -m codex_harness.server
```

## Build the optional container

```powershell
docker compose -f codex_harness/compose.yaml --profile container build
```

The container retains the same stdio-only boundary and publishes no ports.
The repository Docker ignore rules exclude local virtual environments, solution
archives, account bindings, and credentials from every build context.
