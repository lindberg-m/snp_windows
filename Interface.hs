module Interface (getOptions, Options (..)) where

import Options.Applicative

data Options = Options {
   input  :: Maybe String
 , wSize  :: Int
 , wStep  :: Int
 , nohead :: Bool
} deriving (Show, Eq)

getOptions = execParser $ info (helper <*> options)
                     ( fullDesc <>
                       progDesc ("Calculate the mean on a sliding window on a three-column, " ++
                                 "tabseparated file, where col1=chromosome, col2=position and " ++
                                 "col3=data") <>
                       header "snp_windows, a sliding-window calculator for SNP-data")

options :: Parser Options
options = Options <$>
  (optional $ strOption 
     (short 'i' <>
      long "input" <> 
      metavar "PATH" <>
      help "if this opt is absent, input is read from stdin"))
  <*> option auto 
        ( short 's' <>
          long "size" <>
          metavar "INT" <>
          help "set the size of the windows used"
        )
  <*> option auto
        ( short 't' <>
          long "step" <>
          metavar "INT" <>
          help "set the size of steps of windows"
        )
  <*> switch
        ( long "no-header" <>
          help "include the first line, since it's not a header"
        )
  
