# Brain eQTL registry

Each JSON file is one independently selectable QTL dataset. `resolution` separates
bulk tissue from cell-type/subtype evidence and `evidence_family` prevents an
overlapping cohort or a redistributed copy of the same result from being counted
as independent replication. Sample size describes the cited release, not a
guarantee that every gene/cell has that N.

Access is deliberately two-part. A registry entry does not grant access:
`public_summary_statistics` identifies releasable aggregate products, while
`controlled_access` describes protected material. The adapter reports missing or
unapproved inputs as `NOT_RUN_ACCESS_REQUIRED` or `NOT_RUN_INPUT_UNAVAILABLE`;
neither is a PASS.
