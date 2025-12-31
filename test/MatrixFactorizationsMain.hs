module Main (main) where

import Test.Hspec
import MatrixFactorizationsSpec

main :: IO ()
main = hspec $ do
  specQR
  specChol
