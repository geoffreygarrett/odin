#!/bin/bash

# Remove old configuration files
rm -rf ../.clwb/.idea/runConfigurations/*

# List of Bazel commands for which we want to generate configurations
commands=("run" "test" "build" "query" "info" "coverage")

# Default target kind
default_kind="cc_binary"

# Parse blazeInfo.xml and generate configurations
while IFS= read -r line
do
    if [[ $line == *"<entry key=\""* ]]; then
        # Extract target from the line
        target="${line#*<entry key=\"}"
        target="${target%%\" value=\"*}"

        # Determine kind from target name
        kind=$default_kind
        if [[ $target == *":test_"* ]]; then
            kind="cc_test"
        fi

        # Generate a configuration for each command
        for command in "${commands[@]}"
        do
            # Skip running, testing, and coverage of targets that start with "@"
            if [[ $target == @* && $command != "build" && $command != "query" && $command != "info" ]]; then
                continue
            fi

            # Create the configuration file
            cat > ../.clwb/.idea/runConfigurations/${command}__${target//[:\/]/__}.xml << EOF
<component name="ProjectRunConfigurationManager">
  <configuration default="false" name="${command}:${target}" type="BlazeCommandRunConfigurationType" factoryName="Bazel Command" nameIsGenerated="true">
    <blaze-settings handler-id="BlazeCLionRunConfigurationHandlerProvider" kind="$kind" blaze-command="$command">
      <blaze-target>$target</blaze-target>
      <env_state>
        <envs />
      </env_state>
    </blaze-settings>
    <method v="2">
      <option name="Blaze.BeforeRunTask" enabled="true" />
    </method>
  </configuration>
</component>
EOF
        done
    fi
done < ../.clwb/.idea/blazeInfo.xml

## Add this script as a run configuration
#cat > .clwb/.idea/runConfigurations/update_bazel_configs.xml << EOF
#<component name="ProjectRunConfigurationManager">
#  <configuration default="false" name="Update Bazel Configs" type="ShConfigurationType" nameIsGenerated="true">
#    <option name="SCRIPT_TEXT" value="" />
#    <option name="INDEPENDENT_SCRIPT_PATH" value="true" />
#    <option name="SCRIPT_PATH" value="\$PROJECT_DIR\$/update_bazel_configs.sh" />
#    <option name="SCRIPT_OPTIONS" value="" />
#    <option name="INDEPENDENT_SCRIPT_WORKING_DIRECTORY" value="true" />
#    <option name="SCRIPT_WORKING_DIRECTORY" value="\$PROJECT_DIR\$" />
#    <option name="INDEPENDENT_INTERPRETER_PATH" value="true" />
#    <option name="INTERPRETER_PATH" value="/bin/bash" />
#    <option name="INTERPRETER_OPTIONS" value="" />
#    <option name="EXECUTE_IN_TERMINAL" value="true" />
#    <option name="EXECUTE_SCRIPT_FILE" value="true" />
#    <envs />
#    <method v="2" />
#  </configuration>
#</component>
#EOF
