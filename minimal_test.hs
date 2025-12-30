{-# LANGUAGE FlexibleContexts #-}
-- Minimal test to verify syntax and types are correct
module MinimalTest where

import Numeric.LinearAlgebra.Sparse
import Data.Sparse.Common
import Control.Monad.Catch
import Control.Monad.Writer

-- Test that eigsQR type signature is correct
testEigsQR :: (MonadThrow m, MonadWriter [String] m) => SpMatrix Double -> m (SpVector Double)
testEigsQR mat = eigsQR 100 False mat

-- Test that qr type signature is correct
testQR :: (MonadThrow m, MonadWriter [String] m) => SpMatrix Double -> m (SpMatrix Double, SpMatrix Double)
testQR = qr

-- Test that givens type signature is correct  
testGivens :: (MonadThrow m, MonadWriter [String] m) => SpMatrix Double -> IxRow -> IxCol -> m (Maybe (SpMatrix Double))
testGivens = givens

main :: IO ()
main = do
  let mat = sparsifySM $ fromListDenseSM 4 [0,1,1,0, 1,0,0,0, 1,0,0,1, 0,0,1,0] :: SpMatrix Double
  
  -- Test eigsQR
  (result, logs) <- runWriterT $ eigsQR 100 False mat
  putStrLn "eigsQR test passed!"
  putStrLn $ "Eigenvalues dimension: " ++ show (dim result)
  
  -- Test qr
  (q, r) <- runWriterT $ qr mat
  putStrLn "QR test passed!"
  putStrLn $ "Q dimensions: " ++ show (dim q)
  putStrLn $ "R dimensions: " ++ show (dim r)
