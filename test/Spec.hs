module Main (main) where

import Test.Hspec
import qualified LibSpec

main :: IO ()
main = hspec LibSpec.spec
