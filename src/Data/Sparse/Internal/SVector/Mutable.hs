{-# language FlexibleContexts #-}
module Data.Sparse.Internal.SVector.Mutable where

import qualified Data.Vector as V 
import qualified Data.Vector.Mutable as VM
-- import qualified Data.Vector.Generic as VG

-- import Control.Arrow (Arrow, (***), (&&&))
import Control.Monad.Primitive

data SMVector m a = SMV { smvDim :: {-# UNPACK #-} !Int,
                          smvIx :: V.Vector Int,
                          smvVal :: VM.MVector (PrimState m) a }


-- instance Show a => Show (SMVector m a) where
--   show (SMV n ix v) = unwords ["SMV (",show n,"),",show nz,"NZ:",show (V.zip ix v)]
--     where nz = V.length ix

-- instance Functor (SMVector m) where
--   fmap f (SMV n ix v) = SMV n ix (fmap f v)

-- instance Foldable (SMVector m) where
--   foldr f z (SMV _ _ v) = foldr f z v

-- instance Traversable (SMVector m) where
--   traverse f (SMV n ix v) = SMV n ix <$> traverse f v





-- unionWithM g z (SMV n1 ixu u) (SMV n2 ixv v) = do
--   vm <- VM.new n
--   (vm', nfin) <- go ixu u ixv v 0 vm
--   let vm'' = VM.take nfin vm'
--   return vm''
--   where
--     n = min n1 n2
--     go iu u_ iv v_ i vm
--           | (VG.null u_ && VG.null v_) || i == n = return (vm , i)
--           | VG.null u_ = do
--                  VM.write vm i (V.head iv, g z (VG.head v_))
--                  go iu u_ (V.tail iv) (VG.tail v_) (i + 1) vm
--           | VG.null v_ = do
--                  VM.write vm i (V.head iu, g (VG.head u_) z)
--                  go (V.tail iu) (VG.tail u_) iv v_ (i + 1) vm
--           | otherwise =  do
--            let (u, us) = headTailG u_
--                (v, vs) = headTailG v_
--                (iu1, ius) = headTail iu
--                (iv1, ivs) = headTail iv
--            if iu1 == iv1 then do VM.write vm i (iu1, g u v)
--                                  go ius us ivs vs (i + 1) vm
--                          else if iu1 < iv1 then do VM.write vm i (iu1, g u z)
--                                                    go ius us iv v_ (i + 1) vm
--                                            else do VM.write vm i (iv1, g z v)
--                                                    go iu u_ ivs vs (i + 1) vm


-- -- * helpers

-- -- both :: Arrow arr => arr b c -> arr (b, b) (c, c)
-- -- both f = f *** f

-- headTailG = VG.head &&& VG.tail

-- headTail = V.head &&& V.tail
