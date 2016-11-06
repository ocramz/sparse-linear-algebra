module Data.Sparse.Internal.CSR where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VM

{-| Compressed Row Storage specification :

   http://netlib.org/utk/people/JackDongarra/etemplates/node373.html

   The compressed row storage (CRS) format puts the subsequent nonzeros of the matrix
   rows in contiguous memory locations. Assuming we have a nonsymmetric sparse matrix
   $A$, we create three vectors: one for floating point numbers (val) and the other
   two for integers (col_ind, row_ptr).

   The val vector stores the values of the nonzero elements of the matrix $A$ as
   they are traversed in a row-wise fashion.
 
   The col_ind vector stores the column indexes of the elements in the val vector,
   that is, if val(k)=a_{i,j}, then  col_ind(k)=j$.

   The row_ptr vector stores the locations in the val vector that start a row;
   that is, if  val(k)=a_{i,j}, then row_ptr(i) <= k < row_ptr(i+1)

-}

data CsrMatrix a =
  CsrMatrix { csrVal :: VU.Vector a,
              csrColInd :: VU.Vector Int,
              csrRowPtr :: VU.Vector Int,
              csrNnz :: {-# UNPACK #-} !Int,
              csrNrows :: {-# UNPACK #-} !Int,
              csrNcols :: {-# UNPACK #-} !Int } deriving Eq

