module TestWindows where
import Snp_parser
import Snp_windows

import Control.Monad (liftM)
import Data.Either
import Control.Monad.Trans.Except
import qualified Data.Text.Lazy as T
import qualified Data.Text.Lazy.IO as TIO
import Data.List
import Snp_windows 

testDat_Snps = [
    SNP (T.pack "chrom_1") 1 0.424
  , SNP (T.pack "chrom_1") 6 0.123
  , SNP (T.pack "chrom_1") 24 0.13
  , SNP (T.pack "chrom_1") 49 0.1023
  , SNP (T.pack "chrom_1") 104 0.943
  , SNP (T.pack "chrom_1") 599 0.904123
--  , SNP "chrom_1" 100043 0.332
--  , SNP "chrom_1" 1000434 0.123
  , SNP (T.pack "chrom_2") 43 0.23
  , SNP (T.pack "chrom_2") 43 1.0
  , SNP (T.pack "chrom_2") 43 0
  , SNP (T.pack "chrom_2") 43 0.2
  ]

testConfig = WindowConfig 100 10 
testState = initWindows testConfig


test :: FilePath -> IO ()
test p = do
  input <- liftM (rights . (map runExcept) . parseSNPlist . tail . T.lines) $ 
             TIO.readFile p
  let chromWindows = calcWindows (WindowConfig { windowSize = 1000, windowStep = 200 })
                                  input 
      windowStats = concat . map (\(chrom, windows) -> map (meanWindow chrom) windows) 
  mapM_ print $ windowStats chromWindows

type Chrom = T.Text

data WindowStats = WindowStats {
    samples :: Int
  , mean    :: Double
  , start   :: Int
  , end     :: Int
  , chrom   :: Chrom
}

instance Show WindowStats where
  show (WindowStats s m st e c) =
     show c ++ "\t" ++ show st ++ "\t" ++ show e ++ "\t" ++ show s ++"\t"++ show m

meanWindow :: Chrom -> Window Double -> WindowStats
meanWindow c w = WindowStats (round n) (sum / n) (Snp_windows.start w) (Snp_windows.end w) c
  where
    (n, sum) = foldl' (\(n', s) x -> (n' + 1, s + x)) (0.0,0.0) (wData w)
