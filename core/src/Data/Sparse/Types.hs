-----------------------------------------------------------------------------
-- |
-- Module      :  Data.Sparse.Types
-- Copyright   :  (c) Marco Zocca 2017
-- License     :  BSD3 (see the file LICENSE)
--
-- Maintainer  :  zocca marco gmail
-- Stability   :  experimental
-- Portability :  portable
--
-- Typeclasses for linear algebra and related concepts
--
-----------------------------------------------------------------------------
module Data.Sparse.Types where

type Rows = Int
type Cols = Int

type IxRow = Int
type IxCol = Int

-- * Lexicographic ordering types

type LexIx = Int

data LexOrd = RowsFirst | ColsFirst deriving (Eq, Show)
