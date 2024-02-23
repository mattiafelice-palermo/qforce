# A scheduler file should be composed of the following parts
# 1. Scheduler settings
#  - Nodes
#  - number of tasks
#  - ram specification
#  - etc.
# 2. Preparatory phase
#  - module activation
#  - conda activation
#  - handling of scratch dir
#  - etc
#  3. Run command
#  - Can either use default command or user customized commands
#  4. Finalization
#  - handling of scratch dir
#  - cleanup
#  - send signal for any dependent job?
#  - etc
