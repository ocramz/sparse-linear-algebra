module Data.Sparse.Common
       ( module X,
         svToSM, outerProdSV, (><), toSV, extractCol, extractRow,
         extractDiagonalDSM,
         matVec, (#>), vecMat, (<#)) where

import Data.Sparse.Utils as X
import Data.Sparse.SpMatrix as X
import Data.Sparse.SpVector as X

import Numeric.Eps as X
import Numeric.LinearAlgebra.Data as X
import Numeric.LinearAlgebra.Class as X

import Numeric.LinearAlgebra.Sparse.IntMap as X

import qualified Data.IntMap as IM

-- | promote a SV to SM
svToSM :: SpVector a -> SpMatrix a
svToSM (SV n d) = SM (n, 1) $ IM.singleton 0 d


-- ** Outer vector product

outerProdSV, (><) :: Num a => SpVector a -> SpVector a -> SpMatrix a
outerProdSV v1 v2 = fromListSM (m, n) ixy where
  m = dim v1
  n = dim v2
  ixy = [(i,j, x * y) | (i,x) <- toListSV v1 , (j, y) <- toListSV v2]

(><) = outerProdSV



-- |Demote (n x 1) or (1 x n) SpMatrix to SpVector
toSV :: SpMatrix a -> SpVector a
toSV (SM (m,n) im) = SV d $ snd . head $ IM.toList im where
  d | m==1 && n==1 = 1
    | m==1 && n>1 = n 
    | n==1 && m>1 = m
    | otherwise = error $ "toSV : incompatible dimensions " ++ show (m,n)

-- |Extract jth column, and place into SpVector
extractCol :: SpMatrix a -> IxCol -> SpVector a
extractCol m j = toSV $ extractColSM m j      


-- |Extract ith row, and place into SpVector
extractRow :: SpMatrix a -> IxRow -> SpVector a
extractRow m i = toSV $ extractRowSM m i







-- | Extract the diagonal as a SpVector (with default 0)
extractDiagonalDSM :: Num a => SpMatrix a -> SpVector a
extractDiagonalDSM mm = fromListDenseSV n $ foldr ins [] ll  where
  ll = [0 .. n - 1]
  n = nrows mm
  ins i acc = mm@@(i,i) : acc




-- ** Matrix action on a vector

{- 
FIXME : matVec is more general than SpVector's :

\m v -> fmap (`dot` v) m
  :: (Normed f1, Num b, Functor f) => f (f1 b) -> f1 b -> f b
-}



-- |Matrix-on-vector
matVec, (#>) :: Num a => SpMatrix a -> SpVector a -> SpVector a
matVec (SM (nr, nc) mdata) (SV n sv)
  | nc == n = SV nr $ fmap (`dot` sv) mdata
  | otherwise = error $ "matVec : mismatching dimensions " ++ show (nc, n)

(#>) = matVec

-- |Vector-on-matrix (FIXME : transposes matrix: more costly than `matVec`, I think)
vecMat, (<#) :: Num a => SpVector a -> SpMatrix a -> SpVector a  
vecMat (SV n sv) (SM (nr, nc) mdata)
  | n == nr = SV nc $ fmap (`dot` sv) (transposeIM2 mdata)
  | otherwise = error $ "vecMat : mismatching dimensions " ++ show (n, nr)

(<#) = vecMat  
