#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <program> [args...]"
    exit 1
fi

LOGFILE="memory_log_$(date +%Y%m%d_%H%M%S).dat"
CGROUP_MEM_FILE="/sys/fs/cgroup/memory.current"

if [ ! -f "$CGROUP_MEM_FILE" ]; then
    echo "Error: $CGROUP_MEM_FILE not found"
    exit 1
fi

# Start the program in the background
time "$@" &
PROGRAM_PID=$!

echo "Logging memory to $LOGFILE while running: $@"

# Gnuplot-compatible header
echo "# Time(s)  Memory(MB)" > "$LOGFILE"
START_TIME=$(date +%s)

# Log every second until the program exits
while kill -0 "$PROGRAM_PID" 2> /dev/null; do
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))
    MEM_BYTES=$(cat "$CGROUP_MEM_FILE")
    MEM_MB=$((MEM_BYTES / 1024 / 1024))
    echo "$ELAPSED $MEM_MB" >> "$LOGFILE"
    sleep 1
done

echo "Done logging. Output saved to $LOGFILE"
cp "$LOGFILE" "last.dat"