# Test design

1. When running a specific test, use `pixi run cargo test -- test_name`, where `test_name` is the name of the test function to run. This allows for quick iteration on a specific test without running the full test suite.
2. To test for changes more broadly (but still relatively quickly), use `pixi run cargo test`
3. To do the final expensive tests, use `CHECKM2DB=/work/microbiome/db/CheckM2_database/uniref100.KO.1.dmnd pixi run cargo test -- --ignored`
4. Add tests to the appropriate test files in the `test` directory. For example for changes to the analyse subcommand, add tests to `test/test_analyse.rs`. For changes to specific functions, add tests to the end of the file where the function is defined. For example, for changes to `src/skani.rs`, add tests to the end of that file.
5. Update the documentation in docs/ to include any relevant information, including usage instructions and examples.
