# Antenna security boundary

DynSysLab uses the shared `the-antenna` repository through a pinned external
dependency. The tracked integration contains only source code, policy
templates, and placeholder configuration.

Concrete provider bindings, access-key files, credentials, account metadata,
and storage or compute identifiers are local-only. Create local operational
configuration by copying the relevant `config/*.template.json` file, then
replace its placeholders. Those local files are ignored by Git and must not
be pasted into issue trackers, chat, logs, or documentation.

The Antenna policy denies capacity-changing or potentially billable operations
until the user records explicit authorization locally. When capacity limits
are encountered, Antenna should report a practical non-billable workaround or
request authorization; it must not upgrade services on its own.
