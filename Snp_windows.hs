module Snp_windows (calcWindows, Window (..), initWindows, WindowConfig (..), SNP (..) ) where

import Control.Monad.State
import Data.Text.Lazy (Text)
import Data.List (groupBy)

type Pos   = Int
type Chrom = Text

data Window a = Window {
    start :: Pos
  , end   :: Pos
  , wData :: [a]
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


calcWindows :: WindowConfig -> [SNP a] -> [(Chrom, [Window a])]
calcWindows wc snps = map toWindows chrSnps 
  where
    chrSnps         = groupBy (\a b -> chrom a == chrom b) snps
    toWindows snps' = (chrom $ head snps', align (initWindows wc) snps')

{- Populate the window-state with snp-data and return
   a list of populated windows -}
align :: WindowState a -> [SNP a] -> [Window a]
align ws [] = active ws   -- Flush out all remaining active windows and stop
align ws (snp:snps) = doneWindows ++ (align newWindows snps)
  where
    (doneWindows, newWindows) = runState (alignSnp snp) ws

{- Compare the current SNP to the window state:
 Active and inactive windows which has endpoints
 above the current SNP should be returned, active
 windows for which the SNP is overlapping should
 be updated to also carry the current SNP-data, 
 and any inactive windows that have overlapping 
 positions should be updated and moved to 
 the 'active' position.  -}
alignSnp :: SNP a -> State (WindowState a) [Window a]
alignSnp snp = do
  ws <- get
  let
    (passed, remaining)     = span downstreamSnp (active ws)
    (passed', remaining')   = span downstreamSnp (inactive ws)
    (overlaps, remaining'') = span overlapSnp remaining'
    newActive               = map updateWindow (remaining ++ overlaps)
    
    downstreamSnp wndw  = (end wndw) <= (pos snp)     --snp is downstream of window 
    overlapSnp wndw     = ((start wndw) <= (pos snp)) && ((end wndw) >= (pos snp))
    updateWindow wndw   = Window (start wndw) (end wndw) ((snpDat snp) : wData wndw)
  
  put $ WindowState newActive remaining''
  return $ passed ++ passed'

{- Generate a window state which has an empty list as -
 - active windows, and an infinite list of inactive   -
 - windows with sizes based on the passed config      -}
initWindows :: WindowConfig -> WindowState a
initWindows cfg = WindowState { active = [], inactive = wInactive }
  where
    wInactive = zipWith (\a b -> Window a b [])
      (map (+1) [0, step ..])
      [size, (step + size) .. ]
    size = windowSize cfg
    step = windowStep cfg
