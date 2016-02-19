import Interface
import Snp_parser
import Calc_windows

import System.IO                   (hGetContents, stdin)
import System.Exit                 (die)
import Control.Monad               (liftM)
import Control.Monad.Trans.Except  (runExcept)
import Data.Either                 (rights, isRight)
import Data.List                   (foldl')
import qualified Data.Text.Lazy as T
import qualified Data.Text.Lazy.IO as TIO

main :: IO ()
main = do 
  opts <- getOptions 
  let getContent = case input opts of
         Nothing  -> TIO.hGetContents stdin
         (Just x) -> TIO.readFile x
      toLines = if nohead opts
                then T.lines
                else tail . T.lines
      config = WindowConfig (wSize opts) (wStep opts)

  txt <- liftM toLines getContent
  printResults  $ getSNPs config (parseSNPs txt)

data WindowStats = WindowStats {
    samples :: Int
  , mean    :: Double
  , start   :: Int
  , end     :: Int
  , chrom   :: Chrom
}

instance Show WindowStats where
  show (WindowStats s m st e c) =
     T.unpack c ++ "\t" ++ show st ++ "\t" ++ show e ++ "\t" ++ show s ++"\t"++ show m

-- (list of chromosomes-windows, failed parsings (starting with a left value))
type SNPparseResults a = ([(Chrom, [Window a])], [Either String (SNP a)])

printResults :: SNPparseResults Double -> IO ()
printResults (res, failed) = do
   mapM_ print . concat $ map getRes res
   mapM_ catchError failed
   where
   --  getRes :: (Chrom, [Window Double]) -> [WindowStats]
     getRes (chr, w) = map (meanWindow chr) w
     catchError (Left err) = die $ "\nError: \n " ++ err ++ "\n"
     catchError _          = return ()

{- 
 - Provided a config for window-sizes, generate a tuple
 - with all windows in as first element, and snps
 - not mapped to windows in the second element - starting
 - with the first error code from the snp-parsing
 -}

getSNPs :: WindowConfig -> [ParsedSNP a] -> SNPparseResults a
getSNPs wc = mapFirst (calcWindows wc) . span isRight . map runExcept
  where
    mapFirst f (x,y) = (f $ rights x, y)

meanWindow :: Chrom -> Window Double -> WindowStats
meanWindow c w = WindowStats {
    samples = (round n),
    mean = (sum / n),
    Main.start = (Calc_windows.start w),
    Main.end = (Calc_windows.end w),
    Main.chrom = c}
  where
    (n, sum) = foldl' (\(n', s) x -> (n' + 1, s + x)) (0.0,0.0) (wData w)
