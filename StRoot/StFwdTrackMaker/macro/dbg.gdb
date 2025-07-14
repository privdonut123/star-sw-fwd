
# Set breakpoints (example breakpoints, adjust as needed)
set breakpoint pending on 

break StFwdTrackMaker::Make
commands
  printf ">>>>>>>>>>>>>>> MAKE >>>>>>>>>>>>>>>>>>>>>>>\n"
  call printMEM()
  list
end

break FwdTracker.h:998
commands
  printf ">>>>>>>>>>>>>>> FwdTracker::doTrackFitting >>>>>>>>>>>>>>>>>>>>>>>\n"
  call printMEM()
  list
end

# Start the program
run