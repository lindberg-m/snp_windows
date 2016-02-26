module Calc_windows 
    (calcWindows, 
     calcWindows',
     WindowStats (..),
     Window (..), 
     Chrom,
     WindowConfig (..), 
     SNP (..) ) where

import Control.Monad.Trans         (lift)
import Control.Monad.Trans.State   (get, put, runState, evalState, State)
import Control.Monad.Trans.Reader  (ReaderT, ask, runReaderT, runReader, Reader)
import Data.Text.Lazy              (Text, unpack)
import Data.Monoid
import Data.List                   (groupBy)

type Pos   = Int
type Chrom = Text

data Window a = Window {
    start :: Pos
  , end   :: Pos
  , wData :: WindowStats a
} deriving (Show, Eq)

data WindowConfig = WindowConfig {
    windowSize :: Int
  , windowStep :: Int
} deriving (Show, Eq)

data WindowState a = WindowState {
    active   :: [Window a]
  , inactive :: [Window a]
} deriving (Eq)

instance Show a => Show (WindowState a) where
  show (WindowState act inact) =
    "WindowState: Active {" ++ (show act) ++ "}"

data SNP a = SNP {
    chrom  :: Chrom
 ,  pos    :: Pos
 ,  snpDat :: a
} deriving (Show, Eq)

data WindowStats a = WindowStats {
    windowSamples :: Int
  , windowDat     :: a
} deriving (Show, Eq)

instance Monoid a => Monoid (WindowStats a) where
  mempty  = WindowStats 0 mempty
  (WindowStats x a) `mappend` (WindowStats y b) = WindowStats (x + y) (a <> b)

calcWindows :: Num a => [SNP a] -> Reader WindowConfig [(Chrom, [Window (Sum a)])]
calcWindows snps = runReaderT (calcWindows' snps) Sum

calcWindows' :: Monoid b => [SNP a] -> ReaderT (a -> b) (Reader WindowConfig) [(Chrom, [Window b])]
calcWindows' snps = do
  start_windows <- lift initWindows
  readValue     <- ask
  let chrSnps      = groupBy (\a b -> chrom a == chrom b) snps
      toWindows xs = (chrom $ head xs, align' xs start_windows)
      align' xs    = evalState (runReaderT (align xs) readValue)
  return $ map toWindows chrSnps


align :: Monoid b => [SNP a] -> ReaderT (a -> b) (State (WindowState b)) [Window b]
align [] = lift (get >>= \a -> return $ active a) -- Flush out all remaining active windows and stop
align (snp:snps) = do
  alignedSnp <- alignSnp snp
  newAligned <- align snps
  return $ alignedSnp ++ newAligned

{- Compare the current SNP to the window state:
 Active and inactive windows which has endpoints
 above the current SNP should be returned, active
 windows for which the SNP is overlapping should
 be updated to also carry the current SNP-data, 
 and any inactive windows that have overlapping 
 positions should be updated and moved to 
 the 'active' position.  -}
alignSnp :: Monoid b => SNP a -> ReaderT (a -> b) (State (WindowState b)) [Window b]
alignSnp snp = do
  adder <- addToWindow
  ws    <- lift get
  let
    (passed, remaining)     = span downstreamSnp (active ws)
    (passed', remaining')   = span downstreamSnp (inactive ws)
    (overlaps, remaining'') = span overlapSnp remaining'
    newActive               = map updateWindow (remaining ++ overlaps)
    
    downstreamSnp wndw  = (end wndw) <= (pos snp)     --snp is downstream of window 
    overlapSnp wndw     = ((start wndw) <= (pos snp)) && ((end wndw) >= (pos snp))
    updateWindow wndw   = adder (snpDat snp) wndw
  
  lift . put $ WindowState newActive remaining''
  return $ passed ++ passed'

addToWindow :: (Monoid b, Monad m) => ReaderT (a -> b) m (a -> Window b -> Window b)
addToWindow = do
  f <- ask
  return $ \x wdw -> Window (start wdw) (end wdw) (WindowStats (1 + (windowSamples $ wData wdw)) 
                     ((f x) <> (windowDat $ wData wdw)))

{- Generate a window state which has an empty list as -
 - active windows, and an infinite list of inactive   -
 - windows with sizes based on the passed config      -}
initWindows :: (Monoid a, Monad m) => ReaderT WindowConfig m (WindowState a)
initWindows = do
  cfg <- ask
  return $ let
    wInactive = zipWith (\a b -> Window a b mempty)
      (map (+1) [0, step ..])
      [size, (step + size) ..]
    size = windowSize cfg
    step = windowStep cfg
    in WindowState { active = [], inactive = wInactive }

