
<details>
<summary>Incremental Sync output</summary>

```
Updating VCS...
Running Bazel info...
Command: /usr/bin/bazel info --tool_tag=ijwb:CLion --curses=no --color=yes --progress_in_terminal_title=no --
Computing VCS working set...
12:41:29	Initial directory update started
12:41:29	Initial directory update finished
12:41:29	Sync started
12:41:32	Sync finished with errors
```

</details>

<details>
<summary>Sub-task Sync readout (relevant redacted version)</summary>

```markdown
Syncing project: Sync (incremental)...
Updating VCS...
Running Bazel info...
Command: /usr/bin/bazel info --tool_tag=ijwb:CLion --curses=no --color=yes --progress_in_terminal_title=no --
Computing VCS working set...

Command: /tmp/xcode_cquery1099.sh
Loading: 
Loading: 
Loading: 0 packages loaded
Analyzing: target @bazel_tools//tools/osx:current_xcode_config (1 packages loaded, 0 targets configured)
INFO: Analyzed target @bazel_tools//tools/osx:current_xcode_config (2 packages loaded, 3 targets configured).
INFO: Found 1 target...

Command: /usr/bin/bazel run --tool_tag=ijwb:CLion -- @bazel_tools//tools/osx:xcode-locator None
Loading: 
Loading: 
Loading: 0 packages loaded
Analyzing: target @bazel_tools//tools/osx:xcode-locator (0 packages loaded, 0 targets configured)
INFO: Analyzed target @bazel_tools//tools/osx:xcode-locator (2 packages loaded, 9 targets configured).
INFO: Found 1 target...
INFO: Running command line: bazel-bin/external/bazel_tools/tools/osx/xcode-locator None
xcode_locator should not be invoked on non-darwin systems

Sync failed
```

</details>

<details>
<summary>Sub-task Sync readout (long version for transparency)</summary>

```markdown
Syncing project: Sync (incremental)...
Updating VCS...
Running Bazel info...
Command: /usr/bin/bazel info --tool_tag=ijwb:CLion --curses=no --color=yes --progress_in_terminal_title=no --
Computing VCS working set...
Sync targets from project view targets:
//tests/...
//examples/...
//src/...

Running Bazel info...
Command: /usr/bin/bazel info --tool_tag=ijwb:CLion --curses=no --color=yes --progress_in_terminal_title=no --
Building 3 Bazel targets...
Build invocation result: success
Parsing build outputs...
Total rules: 40, new/changed: 0, removed: 0
Reading IDE info result...
Reading IDE info result...
Updating target map...
Loaded 0 aspect files, total size 0kB
Prefetching genfiles...
Target map size: 37
Prefetching output artifacts...
Prefetching files...
Refreshing files...
Computing directory structure...
Initializing project SDKs...
Committing project structure...
Workspace has 0 libraries
Workspace has 2 modules
Updating in-memory state...
Command: /tmp/xcode_cquery1099.sh
Loading:
Loading:
Loading: 0 packages loaded
Analyzing: target @bazel_tools//tools/osx:current_xcode_config (1 packages loaded, 0 targets configured)
INFO: Analyzed target @bazel_tools//tools/osx:current_xcode_config (2 packages loaded, 3 targets configured).
INFO: Found 1 target...
INFO: Elapsed time: 0.300s, Critical Path: 0.00s
INFO: 0 processes.
INFO: Build completed successfully, 0 total actions
Command: /usr/bin/bazel run --tool_tag=ijwb:CLion -- @bazel_tools//tools/osx:xcode-locator None
Loading:
Loading:
Loading: 0 packages loaded
Analyzing: target @bazel_tools//tools/osx:xcode-locator (0 packages loaded, 0 targets configured)
INFO: Analyzed target @bazel_tools//tools/osx:xcode-locator (2 packages loaded, 9 targets configured).
INFO: Found 1 target...
[0 / 2] [Prepa] BazelWorkspaceStatusAction stable-status.txt
Target @bazel_tools//tools/osx:xcode-locator up-to-date:
bazel-bin/external/bazel_tools/tools/osx/xcode-locator
INFO: Elapsed time: 0.186s, Critical Path: 0.01s
INFO: 1 process: 1 internal.
INFO: Build completed successfully, 1 total action
INFO: Running command line: bazel-bin/external/bazel_tools/tools/osx/xcode-locator None
xcode_locator should not be invoked on non-darwin systems
3 unique C configurations, 21 C targets
Sync finished

Timing summary:
BlazeInvocation: 2.3s, Prefetching: 0ms, Other: 538ms
Sync failed
```

</details>
