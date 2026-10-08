# Test design

1. When running a specific test, use `pixi run cargo test -- test_name`, where `test_name` is the name of the test function to run. This allows for quick iteration on a specific test without running the full test suite.
2. To test for changes more broadly (but still relatively quickly), use `pixi run cargo test`
3. To do the final expensive tests, submit them to the PBS queue from the repo root: `bash tests/run_expensive_tests_at_cmr.sh`. To run them directly (e.g. inside an interactive job), use `bash tests/run_expensive_tests.sh`.
4. Add tests to the appropriate test files in the `test` directory. For example for changes to the analyse subcommand, add tests to `test/test_analyse.rs`. For changes to specific functions, add tests to the end of the file where the function is defined. For example, for changes to `src/skani.rs`, add tests to the end of that file.
5. Update the documentation in docs/ to include any relevant information, including usage instructions and examples.
