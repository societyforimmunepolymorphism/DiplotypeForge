# Vendored HFX JSON Schema — provenance

`hfx.schema.json` is the authoritative Haplotype Frequency eXchange (HFX) schema,
vendored **verbatim** from the Society for Immune Polymorphism (SIP) repository so
that DiplotypeForge has **no runtime dependency** on the upstream source.

| Field | Value |
|-------|-------|
| Upstream repository | `societyforimmunepolymorphism/hfx` |
| Upstream URL | https://github.com/societyforimmunepolymorphism/hfx |
| Integrated change | PR #10 (`frequencyDataLocation` correction) |
| PR URL | https://github.com/societyforimmunepolymorphism/hfx/pull/10 |
| Merge commit SHA | `a434a217ff7c8f13a73d04b7635e34e76b3ef310` |
| Document version | `0.1.1` (unchanged by the correction) |
| Canonical SHA-256 | `7609fef72ebf381d36d1bfc49c951d6138717df5fb2ead7ff40a7176c4031bd3` |

The **canonical SHA-256** is computed over the schema's semantic content
(`json.dumps(schema, sort_keys=True, separators=(",", ":"))`), so it is invariant
to pretty-printing and key ordering. Reproduce it with:

```python
from diplotypeforge.adapters.hfx import schema_checksum, schema_provenance
print(schema_checksum())     # 7609fef7...
print(schema_provenance())   # full manifest-ready record
```

This record is also written into every run `manifest.json` under the
`hfx_schema` key.

## The correction integrated by PR #10

Earlier DiplotypeforForge builds applied an in-memory `oneOf` → `anyOf`
relaxation to `frequencyDataLocation` because the previous upstream schema's
`oneOf` branches were not mutually exclusive. In the corrected schema the two
branches are mutually exclusive:

- branch 1 is `{"const": "inline"}`;
- branch 2 is a URI string constrained by `pattern: ^(https?|ftps?|file)://.+`.

The literal `"inline"` matches only branch 1, and any supported URI matches only
branch 2, so `oneOf` is now internally consistent. The runtime workaround has
been **removed**; the schema is validated verbatim with no mutation.

## Updating this schema

If a future authoritative release changes the schema:

1. Fetch the new `hfx.schema.json` from the upstream repository at a specific
   commit.
2. Replace this file.
3. Update `HFX_SCHEMA_UPSTREAM_COMMIT` / `HFX_SCHEMA_VERSION` in
   `src/diplotypeforge/adapters/hfx.py`.
4. Update the SHA-256 and commit here.
5. Run `pytest tests/test_hfx_schema_regression.py`.
