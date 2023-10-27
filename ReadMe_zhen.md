# Notes of working on FttSlowMaker testing branch

Builder: Zhen Wang

## Setup the code
Perform a sparse checkout and hop onto the `fwd-tracking` branch
```sh
git clone --no-checkout git@github.com:jdbrice/star-sw-1.git
cd star-sw-1
git config core.sparseCheckout true

touch .git/info/sparse-checkout
echo "StDb" >> .git/info/sparse-checkout
echo "StarDb" >> .git/info/sparse-checkout
echo "StRoot/RTS" >> .git/info/sparse-checkout
echo "StRoot/StBFChain" >> .git/info/sparse-checkout
echo "StRoot/StEvent" >> .git/info/sparse-checkout
echo "StRoot/StFst*" >> .git/info/sparse-checkout
echo "StRoot/StFtt*" >> .git/info/sparse-checkout
echo "StRoot/StFwd*" >> .git/info/sparse-checkout
echo "StRoot/StFcs*" >> .git/info/sparse-checkout

git checkout ftt-slow-sim
```

##How to run the StFttSimuMaker


