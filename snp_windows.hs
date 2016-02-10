import Control.Monad.Trans
import Control.Monad.Reader
import Control.Monad.State

main :: IO ()
main = undefined


type Pos   = Int
type Chrom = String

data Window a = Window {
    start :: Pos
  , end   :: Pos
  , wData :: [a]
}

data WindowConfig = WindowConfig {
    windowSize :: Int
  , windowStep :: Int
}

data WindowState a = WindowState {
    wsChrom  :: Chrom
  , active   :: [Window a]
  , inactive :: [Window a]
}

data SNP a = SNP {
    chrom  :: Chrom
 ,  pos    :: Pos
 ,  snpDat :: a
}

--  At the time of initiation, the current chromosome is unknown, thus a 
--  window-state is not given (Nothing) to the function alignSnp
--
-- alignSnp :: Reader WindowConfig (SNP a -> State (WindowState a) [Window a])
alignSnp :: WindowConfig -> SNP a -> Maybe (WindowState a) -> ((WindowState a), [Window a])
alignSnp wc snp (Nothing) = let 
  windowState = initWindows (windowSize wc) (windowStep wc) (chrom snp)
  in alignSnp' snp windowState

alignSnp wc snp (Just ws)
  | (chrom snp) == (wsChrom ws) = alignSnp' snp ws
  | otherwise = let
      windowState          = initWindows (windowSize wc) (windowStep wc) (chrom snp)
      ret                  = snd $ flushActiveWindows ws
      (windowState', ret') = alignSnp' snp windowState
      in (windowState', (ret ++ ret'))


--align :: Reader WindowConfig (SNP a -> State (WindowState a) [Window a])
alignSnp' :: SNP a -> WindowState a -> ((WindowState a), [Window a])
alignSnp' snp ws = let
  (passed, remaining)     = break downstreamSnp (active ws)
  (passed', remaining')   = break downstreamSnp (inactive ws)
  (overlaps, remaining'') = span overlapSnp remaining'
  newActive               = map updateWindow (remaining ++ overlaps)
  in (WindowState (wsChrom ws) newActive remaining'', passed ++ passed')
    where
      downstreamSnp wndw  = (end wndw) < (pos snp)     --snp is downstream of window 
      overlapSnp wndw     = ((start wndw) <= (pos snp)) && ((end wndw) >= (pos snp))
      updateWindow wndw   = Window (start wndw) (end wndw) ((snpDat snp) : wData wndw)


--flushActiveWindows :: State (WindowState a) [Window a]
--flushAvtiveWindows =
--   ws <- get
--   put $ WindowState (wsChrom ws) [] (inactive ws)
--   return $ active ws
flushActiveWindows :: WindowState a -> (WindowState a, [Window a])
flushActiveWindows ws = let
  ws' = WindowState (wsChrom ws) [] (inactive ws)
  ret = active ws
  in (ws', ret)

data Overlap = Upstream | Downstream | Overlapping
overlap :: Pos -> Pos -> Pos -> Overlap
overlap x start end
  | x < start = Upstream
  | x > end   = Downstream
  | otherwise = Overlapping


initWindows :: Int -> Int -> Chrom -> WindowState a
initWindows size step x = let 
  wActive   = []
  wInactive = zipWith (\a b -> Window a b [])
    (map (+1) [0, step ..])
    [size, (step + size) .. ]
  in WindowState x wActive wInactive



