module Data.Sparse.Types where

type Rows = Int
type Cols = Int

type IxRow = Int
type IxCol = Int

-- * Lexicographic ordering types

type LexIx = Int

data LexOrd = RowsFirst | ColsFirst deriving (Eq, Show)
