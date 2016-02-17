module TestWindows where

import Snp_parser
import Snp_windows

import Control.Monad (liftM)
import Data.Maybe
import Data.Either
import Control.Monad.Trans.Except
import System.Exit
import System.IO
import qualified Data.Text.Lazy as T
import qualified Data.Text.Lazy.IO as TIO
import Data.List

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
  input <- liftM (rights . (map runExcept) . parseSNPs . tail . T.lines) $ 
             TIO.readFile p
  let chromWindows = calcWindows (WindowConfig { windowSize = 1000, windowStep = 200 })
                                  input 
      windowStats = concat . map (\(chrom, windows) -> map (meanWindow chrom) windows) 
  mapM_ print $ windowStats chromWindows

--data ParseError = UnsortedSNPs 
--                | PositionUnparsed String
--                | ValueUnparsed String
--                | TooFewFields Int Int
--                | ParseError String
--                deriving (Show, Eq)

reportSNPs :: [Except ParseError (SNP a)] -> IO [SNP a]
reportSNPs parsed_snps = do
  m_snps <- sequence . (map reportSNPerror) $ enumerate parsed_snps
  return $ catMaybes m_snps
  where
    enumerate = zip [1..]


-- take line-number and parsed SNP and performs IO to report parse-error to user.
reportSNPerror :: (Int, Except ParseError (SNP a)) -> IO (Maybe (SNP a))
reportSNPerror (i,x) = 
  case runExcept x of
    Left UnsortedSNPs -> do 
      printErr $ "The SNPs are not in order (line " ++ show i ++
                  " and " ++ show (i + 1) ++ ")."
    Left (PositionUnparsed e) -> do
      printErr $ "The position-field could not be parsed at line " ++ show i ++ ": " ++ e
    Left (ValueUnparsed e) -> do
      printErr $ "Could not parse the value-field(s) at line " ++ show i ++": " ++ e
    Left (TooFewFields a b) -> do
      printErr $ "Missing fields at line " ++ show i ++ ". " ++ show a ++ " fields present when " ++
                     show b ++ " fields is the minimum required"
    Left (ParseError e) -> do
      printErr $ "Couldn't parse line " ++ show i ++ ": " ++ e
    Right x -> return $ Just x
        
  where
    printErr e = do 
      hPutStrLn stderr $ "Error: " ++ e
      exitFailure
      return Nothing
    

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
