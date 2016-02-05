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
 ,  snpDat :: [a]
}

-- alignSnp :: Reader WindowConfig (SNP a -> State (WindowState a) [Window a])
alignSnp :: WindowConfig -> SNP a -> Maybe (WindowState a) -> ((WindowState a), [Window a])
alignSnp wc snp (Nothing) = let 
  windowState = initWindows (windowSize wc) (windowStep wc) (chrom snp)
  in align snp windowState

alignSnp wc snp (Just ws)
  | (chrom snp) == (wsChrom ws) = align snp ws
  | otherwise = (initWindows (windowSize wc) (windowStep wc) (chrom snp), 
                 flushActiveWindows ws)


align :: SNP a -> WindowState a -> ((WindowState a), [Window a])
align = undefined


flushActiveWindows :: WindowState a -> [Window a]
flushActiveWindows = undefined

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



