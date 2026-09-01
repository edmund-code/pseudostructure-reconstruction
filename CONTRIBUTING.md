# Contributing

Keep research inputs and generated outputs outside Git. New workflows must accept explicit paths
or documented configuration, never a user-specific directory. Add synthetic coverage for reusable
logic, retain scientific caveats beside result summaries, and run the repository hygiene checks
before committing.

Changes to `analysis/` must preserve the v4 segmentation centroid guard. Changes to
`segmentation/` must retain its CPU unit-test path and must not silently enable expensive model
downloads in CI.
