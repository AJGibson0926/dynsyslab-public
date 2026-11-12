FROM python:3.12-slim

WORKDIR /app
COPY codex_harness /app/codex_harness
RUN pip install --no-cache-dir "mcp>=2,<3"

# MCP stdio is the only enabled transport.  Do not publish a network port.
ENTRYPOINT ["python", "-m", "codex_harness.server"]
