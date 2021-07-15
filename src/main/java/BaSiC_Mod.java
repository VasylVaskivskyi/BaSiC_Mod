/*
MODIFIED VERSION OF BASIC, DECOMPILED FROM BaSiC_.jar
Tingying Peng, Kurt Thorn, Timm Schroeder,
Lichao Wang, Fabian J Theis, Carsten Marr*, Nassir Navab*,
A BaSiC Tool for Background and Shading Correction of Optical Microscopy Images
Nature Communication 8:14836 (2017). doi: 10.1038/ncomms14836.
https://github.com/marrlab/BaSiC
*/
//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by Fernflower decompiler)
//


import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleSingularValueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import cern.jet.math.tdouble.DoublePlusMultSecond;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.NewImage;
import ij.gui.Plot;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;
import java.awt.Font;
import java.awt.TextField;
import java.io.File;
import java.nio.file.Paths;
import java.util.EmptyStackException;


public class BaSiC_Mod implements PlugIn {
    private ImageStack stack;
    private int noOfSlices;
    private int noOfChannels;
    private static String[] shadingEstimationOptions = new String[]{"Skip estimation and use predefined shading profiles", "Estimate shading profiles"};
    private static String[] shadingModelOptions = new String[]{"Estimate flat-field only (ignore dark-field)", "Estimate both flat-field and dark-field"};
    private static String[] parameterSettingOptions = new String[]{"Automatic", "Manual"};
    private static String[] driftOptions = new String[]{"Ignore", "Replace with zero", "Replace with temporal mean"};
    private static String[] correctionOptions = new String[]{"Compute shading and correct images", "Compute shading only"};
    private static String none = "None";

    public BaSiC_Mod() {
    }

    public void run(String arg) {
        this.showDialog();
    }

    void showDialog() {
        GenericDialog gd = new GenericDialog("BaSiC");
        gd.addMessage("BaSiC: A Tool for Background and Shading Correction of Optical Microscopy Images", new Font("SansSerif", 1, 13));
        gd.addMessage("");
        gd.addStringField("Input_stack: ", "", 45);
        gd.addMessage("");
        gd.addRadioButtonGroup("Shading_estimation", shadingEstimationOptions, 1, 2, shadingEstimationOptions[1]);
        gd.addStringField("Flat-field_image_path", "", 45);
        gd.addStringField("Dark-field_image_path", "", 45);
        gd.addMessage("");
        gd.addRadioButtonGroup("Shading_model:", shadingModelOptions, 1, 2, shadingModelOptions[0]);
        gd.addMessage("");
        gd.addRadioButtonGroup("Setting_regularisation_parameters:", parameterSettingOptions, 1, 2, parameterSettingOptions[0]);
        gd.addNumericField("lambda_flat 0-5", 0.5D);
        gd.addNumericField("lambda_dark 0-5", 0.5D);
        gd.addMessage("");
        gd.addRadioButtonGroup("Temporal_drift of baseline (ignore or remove & replace):", driftOptions, 1, 3, driftOptions[0]);
        gd.addMessage("");
        gd.addRadioButtonGroup("Correction_options:", correctionOptions, 1, 2, correctionOptions[0]);
        gd.addMessage("");
        gd.addStringField("Output_dir: ", "", 45);
        gd.addMessage("");
        gd.addMessage("Copyright © 2016 Tingying Peng, Helmholtz Zentrum München and TUM, Germany. All rights reserved.");
        //TextField flat_lambda = (TextField)gd.getNumericFields().get(0);
        //TextField dark_lambda = (TextField)gd.getNumericFields().get(1);
        //dark_lambda.setEnabled(false);
        //flat_lambda.setEnabled(false);
        //gd.addDialogListener(this);
        gd.showDialog();
        if (!gd.wasCanceled()) {
            this.processDialog(gd);
        }
    }

    void processDialog(GenericDialog gd) {
        String stackPath = gd.getNextString();
        String flatPath = gd.getNextString();
        String darkPath = gd.getNextString();
        String outDir = gd.getNextString();
        ImagePlus imp = IJ.openImage(stackPath);
        ImagePlus imp_flat = null;
        ImagePlus imp_dark = null;

        String myShadingEstimationChoice = gd.getNextRadioButton();
        String myShadingModelChoice = gd.getNextRadioButton();
        String myParameterChoice = gd.getNextRadioButton();
        String myDriftChoice = gd.getNextRadioButton();
        String myCorrectionChoice = gd.getNextRadioButton();
        double lambda_flat = gd.getNextNumber();
        double lambda_dark = gd.getNextNumber();

        if (imp.getBitDepth() == 24) {
            IJ.error("Please decompose RGB images into single channels.");
        } else {
            this.noOfSlices = imp.getNSlices();
            if (this.noOfSlices == 1)
            {
                IJ.error("Input must be an image stack, not a single image");
            }
        }

        if (myShadingEstimationChoice.equals(shadingEstimationOptions[0]))
        {
            imp_flat = IJ.openImage(flatPath);
            imp_dark = IJ.openImage(darkPath);
            if (imp_flat.getNSlices() != 1) {
                IJ.error("Flat-field input must be a single image.");
                return;
            }
            if (imp_dark.getNSlices() != 1) {
                IJ.error("Dark-field input must be a single image.");
                return;
            }
        }

        this.exec(imp, imp_flat, imp_dark,
                myShadingEstimationChoice, myShadingModelChoice, myParameterChoice,
                lambda_flat, lambda_dark,
                myDriftChoice, myCorrectionChoice,
                outDir);

    }

    public void exec(ImagePlus imp, ImagePlus imp_flat, ImagePlus imp_dark,
                     String myShadingEstimationChoice, String myShadingModelChoice, String myParameterChoice,
                     double lambda_flat, double lambda_dark,
                     String myDriftChoice, String myCorrectionChoice,
                     String outDir)
    {
        this.stack = imp.getStack();
        int outputWidth = this.stack.getWidth();
        int outputHeight = this.stack.getHeight();
        String title = imp.getTitle();
        BaSiC_Mod.Options myOptions = new BaSiC_Mod.Options();
        if (myShadingEstimationChoice.equals(shadingEstimationOptions[1])) {
            myOptions.shadingEst = true;
        } else if (myShadingEstimationChoice.equals(shadingEstimationOptions[0])) {
            myOptions.shadingEst = false;
        } else {
            IJ.showMessage("Your inputs for Shading_estimation does not fit any of the available options; use default setting instead");
            myOptions.shadingEst = true;
        }

        if (myShadingModelChoice.equals(shadingModelOptions[0])) {
            myOptions.darkfieldEst = false;
        } else if (myShadingModelChoice.equals(shadingModelOptions[1])) {
            myOptions.darkfieldEst = true;
        } else {
            myOptions.darkfieldEst = false;
            IJ.showMessage("Your inputs for Shading_model does not fit any of the available options; use default setting instead");
        }

        if (myParameterChoice.equals(parameterSettingOptions[0])) {
            myOptions.lambdadark_auto = true;
            myOptions.lambda_auto = true;
        } else if (myParameterChoice.equals(parameterSettingOptions[1])) {
            myOptions.lambdadark_auto = false;
            myOptions.lambda_auto = false;
            myOptions.lambda_dark = lambda_dark;
            myOptions.lambda = lambda_flat;
        } else {
            myOptions.lambdadark_auto = true;
            myOptions.lambda_auto = true;
            IJ.showMessage("Your inputs for Setting regularisation parameters does not fit any of the available options; use default setting instead");
        }

        if (myDriftChoice.equals(driftOptions[0])) {
            myOptions.driftOpt = 0;
        } else if (myDriftChoice.equals(driftOptions[1])) {
            myOptions.driftOpt = 1;
        } else if (myDriftChoice.equals(driftOptions[2])) {
            myOptions.driftOpt = 2;
        } else {
            myOptions.driftOpt = 0;
            IJ.showMessage("Your inputs for Temporal drift does not fit any of the available options; use default setting instead");
        }

        if (myCorrectionChoice.equals(correctionOptions[0])) {
            myOptions.imageCorr = true;
        } else if (myCorrectionChoice.equals(correctionOptions[1])) {
            myOptions.imageCorr = false;
        } else {
            myOptions.imageCorr = true;
            IJ.showMessage("Your inputs for Correction options does not fit any of the available options; use default setting instead");
        }

        ImageStack processingStack = new ImageStack(128, 128);

        ImageProcessor ip_flatfield;
        ImageProcessor ip_darkfield;
        for(int j = 1; j <= this.noOfSlices; ++j) {
            ip_flatfield = this.stack.getProcessor(j).convertToFloat();
            ip_flatfield.setInterpolate(true);
            ip_flatfield.setInterpolationMethod(1);
            ip_darkfield = ip_flatfield.resize(128, 128, true);
            processingStack.addSlice(ip_darkfield);
        }

        BaSiC_Mod.Shading processingShading = new BaSiC_Mod.Shading();
        ip_flatfield = null;
        ip_darkfield = null;
        ImageProcessor imageResized;
        if (!myOptions.shadingEst) {
            if (imp_flat != null) {
                imageResized = imp_flat.getProcessor().convertToFloat();
                imageResized.setInterpolate(true);
                imageResized.setInterpolationMethod(1);
                DoubleMatrix2D flatfieldMatrix = this.imageToMatrix(imageResized.resize(128, 128, true));
                flatfieldMatrix.normalize();
                flatfieldMatrix.assign(DoubleFunctions.mult(16384.0D));
                processingShading.flatfield = new ImagePlus("flatfield", this.matrixToImage(flatfieldMatrix));
                ip_flatfield = processingShading.flatfield.getProcessor();
            } else {
                processingShading.flatfield = NewImage.createFloatImage("flatfield", 128, 128, 1, 2);
                ip_flatfield = processingShading.flatfield.getProcessor();
                ip_flatfield.set(1.0D);
            }

            if (imp_dark != null) {
                imageResized = imp_dark.getProcessor().convertToFloat();
                imageResized.setInterpolate(true);
                imageResized.setInterpolationMethod(1);
                processingShading.darkfield = new ImagePlus("darkfield", imageResized.resize(128, 128, true));
                ip_darkfield = processingShading.darkfield.getProcessor();
            } else {
                processingShading.darkfield = NewImage.createFloatImage("darkfield", 128, 128, 1, 2);
                ip_darkfield = processingShading.darkfield.getProcessor();
                ip_darkfield.set(0.0D);
            }
        } else {
            processingShading = this.ShadingCorrection(processingStack, myOptions);
            ip_flatfield = processingShading.flatfield.getProcessor();
            ip_darkfield = processingShading.darkfield.getProcessor();
        }

        ip_flatfield.setInterpolate(true);
        ip_flatfield.setInterpolationMethod(2);
        imageResized = ip_flatfield.resize(outputWidth, outputHeight);
        ip_darkfield.setInterpolate(true);
        ip_darkfield.setInterpolationMethod(2);
        ImageProcessor ip_outputdarkfield = ip_darkfield.resize(outputWidth, outputHeight);
        BaSiC_Mod.Baseline processingBaseline = this.BaselineCorrection(processingStack, myOptions, ip_flatfield, ip_darkfield);
        float[] basefluor = new float[this.noOfSlices];
        float meanbasefluor = 0.0F;

        for(int i = 0; i < this.noOfSlices; ++i) {
            basefluor[i] = (float)processingBaseline.basefluor[i];
            meanbasefluor += basefluor[i];
        }

        meanbasefluor /= (float)this.noOfSlices;
        BaSiC_Mod.Shading outputShading = new BaSiC_Mod.Shading();
        String newTitleFlat = "flatfield_".concat(title);
        outputShading.flatfield = new ImagePlus(newTitleFlat, imageResized);
        SaveImage(outDir, "flatfield", newTitleFlat, outputShading.flatfield);
        if (myOptions.darkfieldEst) {
            String newTitleDark = "darkfield_".concat(title);
            outputShading.darkfield = new ImagePlus(newTitleDark, ip_outputdarkfield);
            SaveImage(outDir, "darkfield", newTitleDark, outputShading.darkfield);
        }

        if (myOptions.driftOpt != 0) {
            IJ.log("Temporal components:");
            float[] x = new float[this.noOfSlices];
            ResultsTable temporal_rt = new ResultsTable();
            Analyzer.setResultsTable(temporal_rt);

            for(int j = 1; j <= this.noOfSlices; ++j) {
                x[j - 1] = (float)j;
                temporal_rt.incrementCounter();
                temporal_rt.addValue("FrameNo", (double)j);
                temporal_rt.addValue("Basefluor", (double)basefluor[j - 1]);
            }

            temporal_rt.show("Temporal components");
            Plot basefluor_plot = new Plot("Basefluor", "FrameNo", "Basefluor", x, basefluor);
            basefluor_plot.show();
        }

        if (myOptions.imageCorr) {
            ImageStack corrected_stack = new ImageStack(imp.getWidth(), imp.getHeight());

            for(int j = 1; j <= this.noOfSlices; ++j) {
                ImageProcessor imageCorrected = imp.getStack().getProcessor(j).duplicate().convertToFloat();
                float[] imageCorrectedPixels = (float[])imageCorrected.getPixels();
                int i;
                if (myOptions.darkfieldEst) {
                    for(i = 0; i < imageCorrected.getWidth() * imageCorrected.getHeight(); ++i) {
                        imageCorrectedPixels[i] = (imageCorrectedPixels[i] - outputShading.darkfield.getProcessor().getf(i)) / outputShading.flatfield.getProcessor().getf(i);
                    }
                } else {
                    for(i = 0; i < imageCorrected.getWidth() * imageCorrected.getHeight(); ++i) {
                        imageCorrectedPixels[i] /= outputShading.flatfield.getProcessor().getf(i);
                    }
                }

                if (myOptions.driftOpt != 0) {
                    if (myOptions.driftOpt == 1) {
                        for(i = 0; i < imageCorrected.getWidth() * imageCorrected.getHeight(); ++i) {
                            imageCorrectedPixels[i] -= basefluor[j - 1];
                        }
                    } else if (myOptions.driftOpt == 2) {
                        for(i = 0; i < imageCorrected.getWidth() * imageCorrected.getHeight(); ++i) {
                            imageCorrectedPixels[i] = imageCorrectedPixels[i] - basefluor[j - 1] + meanbasefluor;
                        }
                    }
                }

                if (imp.getBitDepth() == 8) {
                    imageCorrected = imageCorrected.convertToByte(false);
                } else if (imp.getBitDepth() == 16) {
                    imageCorrected = imageCorrected.convertToShort(false);
                }

                corrected_stack.addSlice(imp.getStack().getShortSliceLabel(j), imageCorrected);
            }

            String newTitleCor = "corrected_".concat(title);
            ImagePlus corrected_imp = new ImagePlus(newTitleCor, corrected_stack);
            SaveImage(outDir, "corrected", newTitleCor, corrected_imp);
        }

    }

    private void SaveImage(String outDir, String imageType, String title, ImagePlus image){
        String fullOutDir = Paths.get(outDir, imageType).toString();
        String outPath = Paths.get(fullOutDir, title).toString();
        boolean res = new File(fullOutDir).mkdirs();
        IJ.saveAs(image, "Tiff", outPath);
    }

    private final BaSiC_Mod.Shading ShadingCorrection(ImageStack stack, BaSiC_Mod.Options myOptions) {
        int stackWidth = stack.getWidth();
        int stackHeight = stack.getHeight();
        int stackDim = stackWidth * stackHeight;
        int stackZSlice = stack.getSize();
        DoubleMatrix2D D = this.stackToMatrix(stack);
        DenseDoubleMatrix2D weight;
        if (myOptions.lambda_auto || myOptions.lambdadark_auto) {
            DoubleMatrix1D meanD = this.matrixMean(D);
            meanD.assign(DoubleFunctions.div(meanD.zSum() / (double)meanD.size()));
            weight = new DenseDoubleMatrix2D(128, 128);
            weight.assign(meanD.reshape(128, 128));
            weight.dct2(true);
            Double abs_W = weight.assign(DoubleFunctions.abs).zSum();
            if (myOptions.lambda_auto) {
                myOptions.lambda = abs_W / 400.0D * 0.5D;
            }

            if (myOptions.lambdadark_auto) {
                myOptions.lambda_dark = abs_W / 400.0D * 0.2D;
            }

            IJ.log("Smooth regularisation parameters:\n");
            if (myOptions.darkfieldEst) {
                IJ.log("λ_flat = " + myOptions.lambda + "; λ_dark = " + myOptions.lambda_dark);
            } else {
                IJ.log("λ_flat = " + myOptions.lambda);
            }
        }

        IJ.log("Compute spatial shading profiles...");
        DoubleMatrix2D Dsorted = new DenseDoubleMatrix2D(stackDim, stackZSlice);

        for(int i = 0; i < D.rows(); ++i) {
            Dsorted.viewRow(i).assign(D.viewRow(i).viewSorted());
        }

        weight = new DenseDoubleMatrix2D(Dsorted.rows(), Dsorted.columns());
        weight.assign(1.0D);
        BaSiC_Mod.DecomposedMatrix D_decompose = new BaSiC_Mod.DecomposedMatrix();

        for(int iter = 1; iter <= 5; ++iter) {
            IJ.log("Reweighting Iteration:" + iter);
            D_decompose = this.inexactAlmRspcaL1(Dsorted, weight, myOptions);
            DoubleMatrix1D XARowMean = this.matrixMean(D_decompose.LowrankComponent.viewDice());
            weight.assign(D_decompose.SparseComponent);

            for(int u = 0; u < weight.columns(); ++u) {
                weight.viewColumn(u).assign(DoubleFunctions.div(XARowMean.getQuick(u) + 1.0E-6D));
            }

            weight.assign(DoubleFunctions.chain(DoubleFunctions.inv, DoubleFunctions.chain(DoubleFunctions.plus(0.1D), DoubleFunctions.abs)));
            weight.normalize();
            weight.assign(DoubleFunctions.mult((double)weight.size()));
        }

        DoubleMatrix1D temp_XA = this.matrixMean(D_decompose.LowrankComponent);
        temp_XA.assign(D_decompose.Offset, DoubleFunctions.minus);
        DoubleMatrix2D flatfieldMatrix = new DenseDoubleMatrix2D(128, 128);
        flatfieldMatrix.assign(temp_XA.reshape(128, 128));
        flatfieldMatrix.normalize();
        flatfieldMatrix.assign(DoubleFunctions.mult(16384.0D));
        DoubleMatrix2D darkfieldMatrix = new DenseDoubleMatrix2D(128, 128);
        darkfieldMatrix.assign(D_decompose.Offset.reshape(128, 128));
        BaSiC_Mod.Shading shadingExp = new BaSiC_Mod.Shading();
        shadingExp.darkfield = new ImagePlus("darkfield", this.matrixToImage(darkfieldMatrix));
        shadingExp.flatfield = new ImagePlus("flatfield", this.matrixToImage(flatfieldMatrix));
        return shadingExp;
    }

    private final BaSiC_Mod.Baseline BaselineCorrection(ImageStack stack, BaSiC_Mod.Options myOptions, ImageProcessor flatfield, ImageProcessor darkfield) {
        int stackWidth = stack.getWidth();
        int stackHeight = stack.getHeight();
        int var10000 = stackWidth * stackHeight;
        int stackZSlice = stack.getSize();
        DoubleMatrix2D D = this.stackToMatrix(stack);
        DoubleMatrix2D flatfieldMatrix = this.imageToMatrix(flatfield);
        flatfieldMatrix.normalize();
        flatfieldMatrix.assign(DoubleFunctions.mult(16384.0D));
        DoubleMatrix2D darkfieldMatrix = this.imageToMatrix(darkfield);
        DoubleMatrix1D baseFluor = new DenseDoubleMatrix1D(D.rows());
        DoubleMatrix2D weight = new DenseDoubleMatrix2D(D.rows(), D.columns());
        weight.assign(1.0D);
        new BaSiC_Mod.DecomposedMatrix();
        if (myOptions.driftOpt != 0) {
            IJ.log("Compute temporal components...");

            for(int iter = 1; iter <= 2; ++iter) {
                IJ.log("Reweighting Iteration:" + iter);
                BaSiC_Mod.DecomposedMatrix D_decompose = this.baseFluorEst(D, weight, flatfieldMatrix, darkfieldMatrix, myOptions);
                baseFluor = this.matrixMean(D_decompose.LowrankComponent.viewDice());
                ((DoubleMatrix1D)baseFluor).assign(D_decompose.Coeff);
                weight.assign(D_decompose.SparseComponent);

                for(int u = 0; u < weight.columns(); ++u) {
                    weight.viewColumn(u).assign(DoubleFunctions.div(((DoubleMatrix1D)baseFluor).getQuick(u) + 1.0E-8D));
                }

                weight.assign(DoubleFunctions.chain(DoubleFunctions.inv, DoubleFunctions.chain(DoubleFunctions.plus(0.1D), DoubleFunctions.abs)));
                weight.normalize();
                weight.assign(DoubleFunctions.mult((double)weight.size()));
            }
        }

        BaSiC_Mod.Baseline baselineExp = new BaSiC_Mod.Baseline();
        baselineExp.basefluor = ((DoubleMatrix1D)baseFluor).toArray();
        return baselineExp;
    }

    private final BaSiC_Mod.DecomposedMatrix inexactAlmRspcaL1(DoubleMatrix2D D, DoubleMatrix2D weight, BaSiC_Mod.Options myOptions) {
        DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(D, false, false);
        double norm_two = svd.norm2();
        DenseDoubleAlgebra Algebra = new DenseDoubleAlgebra();
        double mu = 12.5D / norm_two;
        double mu_bar = mu * 1.0E7D;
        double rho = 1.5D;
        double normF_D = Algebra.normF(D);
        DoubleMatrix2D A_hat = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DoubleMatrix2D E_hat = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DoubleMatrix1D A_offset = new DenseDoubleMatrix1D(D.rows());
        A_offset.assign(0.0D);
        DoubleMatrix1D A_coeff = new DenseDoubleMatrix1D(D.columns());
        A_coeff.assign(1.0D);
        DoubleMatrix2D Y1 = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DoubleMatrix2D Z1 = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DenseDoubleMatrix2D temp_Wmean = new DenseDoubleMatrix2D(128, 128);
        DoubleMatrix2D W_hat = new DenseDoubleMatrix2D(128, 128);
        DenseDoubleMatrix2D W_idct_hat = new DenseDoubleMatrix2D(128, 128);
        double[] D_min = D.getMinLocation();
        double B1_uplimit = D_min[0];
        double B1_offset = 0.0D;
        int iter = 0;
        boolean converged = false;

        while(!converged) {
            ++iter;
            W_idct_hat.assign(W_hat);
            W_idct_hat.idct2(true);
            Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, (DoubleMatrix2D)A_hat);

            int v;
            for(v = 0; v < ((DoubleMatrix2D)A_hat).rows(); ++v) {
                ((DoubleMatrix2D)A_hat).viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
            }

            temp_Wmean.assign(this.matrixMean(this.computeResidual(D, (DoubleMatrix2D)A_hat, E_hat, Y1, mu, 1.0D)).reshape(128, 128));
            temp_Wmean.dct2(true);
            W_hat.assign(temp_Wmean, DoubleFunctions.plus);
            W_hat.assign(this.shrinkageOperator((DoubleMatrix2D)W_hat, myOptions.lambda / (1.0D * mu)));
            W_idct_hat.assign(W_hat);
            W_idct_hat.idct2(true);
            A_hat = Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, (DoubleMatrix2D)A_hat);

            for(v = 0; v < ((DoubleMatrix2D)A_hat).rows(); ++v) {
                ((DoubleMatrix2D)A_hat).viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
            }

            E_hat.assign(this.computeResidual(D, (DoubleMatrix2D)A_hat, E_hat, Y1, mu, 1.0D), DoubleFunctions.plus);
            E_hat.assign(this.shrinkageOperator(E_hat, weight, 1.0D / (1.0D * mu)));
            Z1.assign(D).assign(E_hat, DoubleFunctions.minus);
            A_coeff.assign(this.matrixMean(Z1.viewDice()).assign(DoubleFunctions.div(Z1.zSum() / (double)Z1.size())));
            A_coeff.assign(DoubleFunctions.max(0.0D));
            if (myOptions.darkfieldEst) {
                DoubleMatrix1D B_offset = new DenseDoubleMatrix1D(D.rows());
                DoubleMatrix1D A1_offset = new DenseDoubleMatrix1D(D.rows());
                DoubleMatrix1D A_coeffminus1 = new DenseDoubleMatrix1D(D.columns());
                IntArrayList validAcoeffList = new IntArrayList();
                DoubleArrayList valueList = new DoubleArrayList();
                A_coeffminus1.assign(A_coeff).assign(DoubleFunctions.minus(1.0D));
                A_coeffminus1.getNegativeValues(validAcoeffList, valueList);
                double mean_W_idct = W_idct_hat.zSum() / (double)W_idct_hat.size();
                IntArrayList InMaskList = new IntArrayList();
                IntArrayList OutMaskList = new IntArrayList();
                DoubleMatrix1D compare_W_idct_hat = new DenseDoubleMatrix1D(D.rows());
                compare_W_idct_hat.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.minus(mean_W_idct - 1.0E-6D));
                compare_W_idct_hat.getPositiveValues(InMaskList, valueList);
                compare_W_idct_hat.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.minus(mean_W_idct + 1.0E-6D));
                compare_W_idct_hat.getNegativeValues(OutMaskList, valueList);
                int[] validAcoeff = this.validElements(validAcoeffList);
                int[] InMask = this.validElements(InMaskList);
                int[] OutMask = this.validElements(OutMaskList);
                DoubleMatrix1D B_coeff = new DenseDoubleMatrix1D(validAcoeff.length);
                B_coeff.assign(this.matrixMean(Z1.viewSelection(InMask, validAcoeff).viewDice()));
                B_coeff.assign(this.matrixMean(Z1.viewSelection(OutMask, validAcoeff).viewDice()), DoubleFunctions.minus);
                B_coeff.assign(DoubleFunctions.div(Z1.zSum() / (double)Z1.size()));
                double temp1 = A_coeff.viewSelection(validAcoeff).aggregate(DoubleFunctions.plus, DoubleFunctions.square);
                double temp2 = A_coeff.viewSelection(validAcoeff).zSum();
                double temp3 = B_coeff.zSum();
                double temp4 = A_coeff.viewSelection(validAcoeff).aggregate(B_coeff, DoubleFunctions.plus, DoubleFunctions.mult);
                double temp6 = temp2 * temp3 - (double)validAcoeff.length * temp4;
                if (temp6 != 0.0D) {
                    B1_offset = (temp1 * temp3 - temp2 * temp4) / temp6;
                } else {
                    B1_offset = 0.0D;
                }

                B1_offset = Math.max(B1_offset, 0.0D);
                B1_offset = Math.min(B1_offset, B1_uplimit / Math.max(mean_W_idct, 1.0E-6D));
                B_offset.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.mult(-1.0D * B1_offset)).assign(DoubleFunctions.plus(B1_offset * mean_W_idct));
                double temp5 = A_coeff.viewSelection(validAcoeff).zSum() / (double)validAcoeffList.size();
                A1_offset.assign(W_idct_hat.vectorize()).assign(DoubleFunctions.mult(-1.0D * temp5));
                A1_offset.assign(this.matrixMean(Z1.viewSelection((int[])null, validAcoeff)), DoubleFunctions.plus);
                A_offset.assign(A1_offset).assign(DoubleFunctions.minus(A1_offset.zSum() / (double)A1_offset.size())).assign(B_offset, DoubleFunctions.minus);
                DenseDoubleMatrix2D W_offset = new DenseDoubleMatrix2D(128, 128);
                W_offset.assign(A_offset.reshape(128, 128));
                W_offset.dct2(true);
                W_offset.assign(this.shrinkageOperator((DoubleMatrix2D)W_offset, myOptions.lambda_dark / (10.0D * mu)));
                W_offset.idct2(true);
                A_offset.assign(W_offset.vectorize());
                A_offset.assign(this.shrinkageOperator((DoubleMatrix1D)A_offset, myOptions.lambda_dark / (10.0D * mu)));
                A_offset.assign(B_offset, DoubleFunctions.plus);
            }

            Z1.assign(D).assign((DoubleMatrix2D)A_hat, DoubleFunctions.minus).assign(E_hat, DoubleFunctions.minus);
            double normF_Z1 = Algebra.normF(Z1);
            Y1.assign(Z1, DoublePlusMultSecond.plusMult(mu));
            mu = mu * 1.5D < mu_bar ? mu * 1.5D : mu_bar;
            double stopCriterion = normF_Z1 / normF_D;
            IJ.log("Stop Criterion" + iter + "  " + stopCriterion);
            if (stopCriterion < 1.0E-6D) {
                converged = true;
            }

            if (!converged && (double)iter >= 500.0D) {
                IJ.log("Maximum iteration reached");
                converged = true;
            }
        }

        BaSiC_Mod.DecomposedMatrix D_decomposed = new BaSiC_Mod.DecomposedMatrix();
        D_decomposed.LowrankComponent = (DoubleMatrix2D)A_hat;
        D_decomposed.SparseComponent = E_hat;
        D_decomposed.Offset = A_offset;
        D_decomposed.Offset.assign(W_idct_hat.vectorize(), DoublePlusMultSecond.plusMult(B1_offset));
        D_decomposed.Coeff = A_coeff;
        return D_decomposed;
    }

    private final BaSiC_Mod.DecomposedMatrix baseFluorEst(DoubleMatrix2D D, DoubleMatrix2D weight, DoubleMatrix2D flatfield, DoubleMatrix2D darkfield, BaSiC_Mod.Options myOptions) {
        DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(D, false, false);
        double norm_two = svd.norm2();
        DenseDoubleAlgebra Algebra = new DenseDoubleAlgebra();
        double mu = 12.5D / norm_two;
        double mu_bar = mu * 1.0E7D;
        double rho = 1.5D;
        double normF_D = Algebra.normF(D);
        DoubleMatrix2D A_hat = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DoubleMatrix2D E_hat = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DoubleMatrix1D A_offset = darkfield.vectorize();
        DoubleMatrix1D A_coeff = this.matrixMean(D.viewDice());
        DenseDoubleMatrix2D W_idct_hat = new DenseDoubleMatrix2D(128, 128);
        W_idct_hat.assign(flatfield);
        DoubleMatrix2D Y1 = new DenseDoubleMatrix2D(D.rows(), D.columns());
        DoubleMatrix2D Z1 = new DenseDoubleMatrix2D(D.rows(), D.columns());
        int iter = 0;
        boolean converged = false;

        while(!converged) {
            ++iter;
            Algebra.multOuter(W_idct_hat.vectorize(), A_coeff, A_hat);

            for(int v = 0; v < A_hat.rows(); ++v) {
                A_hat.viewRow(v).assign(DoubleFunctions.plus(A_offset.getQuick(v)));
            }

            E_hat.assign(this.computeResidual(D, A_hat, E_hat, Y1, mu, 1.0D), DoubleFunctions.plus);
            E_hat.assign(this.shrinkageOperator(E_hat, weight, 1.0D / (1.0D * mu)));
            Z1.assign(D).assign(E_hat, DoubleFunctions.minus);
            A_coeff.assign(this.matrixMean(Z1.viewDice()));
            A_coeff.assign(DoubleFunctions.minus(A_offset.zSum() / (double)A_offset.size()));
            A_coeff.assign(DoubleFunctions.max(0.0D));
            Z1.assign(D).assign(A_hat, DoubleFunctions.minus).assign(E_hat, DoubleFunctions.minus);
            double normF_Z1 = Algebra.normF(Z1);
            Y1.assign(Z1, DoublePlusMultSecond.plusMult(mu));
            mu = mu * 1.5D < mu_bar ? mu * 1.5D : mu_bar;
            double stopCriterion = normF_Z1 / normF_D;
            IJ.log("Stop Criterion" + iter + "  " + stopCriterion);
            if (stopCriterion < 1.0E-6D) {
                converged = true;
            }

            if (!converged && (double)iter >= 500.0D) {
                IJ.log("Maximum iteration reached");
                converged = true;
            }
        }

        BaSiC_Mod.DecomposedMatrix D_decomposed = new BaSiC_Mod.DecomposedMatrix();
        D_decomposed.LowrankComponent = A_hat;
        D_decomposed.SparseComponent = E_hat;
        D_decomposed.Offset = A_offset;
        D_decomposed.Coeff = A_coeff;
        return D_decomposed;
    }

    public final int[] validElements(IntArrayList array) {
        int[] myElements = new int[array.size()];

        for(int i = 0; i < array.size(); ++i) {
            myElements[i] = array.getQuick(i);
        }

        return myElements;
    }

    public final DoubleMatrix1D matrixMean(DoubleMatrix2D A) {
        DoubleMatrix1D average = new DenseDoubleMatrix1D(A.rows());

        for(int i = 0; i < A.rows(); ++i) {
            average.setQuick(i, A.viewRow(i).zSum() / (double)A.columns());
        }

        return average;
    }

    public final DoubleMatrix1D matrixMedian(DoubleMatrix2D A) {
        DoubleMatrix1D median = new DenseDoubleMatrix1D(A.rows());

        for(int i = 0; i < A.rows(); ++i) {
            median.setQuick(i, A.viewRow(i).viewSorted().getQuick(A.columns() / 2));
        }

        return median;
    }

    private final DoubleMatrix2D computeResidual(DoubleMatrix2D D, DoubleMatrix2D A_hat, DoubleMatrix2D E_hat, DoubleMatrix2D Y1, double mu, double ent) {
        DoubleMatrix2D Residual = new DenseDoubleMatrix2D(D.rows(), D.columns());
        Residual.assign(Y1).assign(DoubleFunctions.div(mu)).assign(D, DoubleFunctions.plus).assign(A_hat, DoubleFunctions.minus).assign(E_hat, DoubleFunctions.minus).assign(DoubleFunctions.div(ent));
        return Residual;
    }

    private final DoubleMatrix2D shrinkageOperator(DoubleMatrix2D a_offset, double epsilon) {
        DoubleMatrix2D A_shrink = new DenseDoubleMatrix2D(a_offset.rows(), a_offset.columns());
        DoubleMatrix2D A_shrinknegative = new DenseDoubleMatrix2D(a_offset.rows(), a_offset.columns());
        A_shrink.assign(a_offset).assign(DoubleFunctions.minus(epsilon)).assign(DoubleFunctions.max(0.0D));
        A_shrinknegative.assign(a_offset).assign(DoubleFunctions.plus(epsilon)).assign(DoubleFunctions.min(0.0D));
        A_shrink.assign(A_shrinknegative, DoubleFunctions.plus);
        return A_shrink;
    }

    private final DoubleMatrix1D shrinkageOperator(DoubleMatrix1D a_offset, double epsilon) {
        DoubleMatrix1D A_shrink = new DenseDoubleMatrix1D((int)a_offset.size());
        DoubleMatrix1D A_shrinknegative = new DenseDoubleMatrix1D((int)a_offset.size());
        A_shrink.assign(a_offset).assign(DoubleFunctions.minus(epsilon)).assign(DoubleFunctions.max(0.0D));
        A_shrinknegative.assign(a_offset).assign(DoubleFunctions.plus(epsilon)).assign(DoubleFunctions.min(0.0D));
        A_shrink.assign(A_shrinknegative, DoubleFunctions.plus);
        return A_shrink;
    }

    private final DoubleMatrix2D shrinkageOperator(DoubleMatrix2D A, DoubleMatrix2D weight, double epsilon) {
        DoubleMatrix2D A_shrink = new DenseDoubleMatrix2D(A.rows(), A.columns());
        DoubleMatrix2D A_shrinknegative = new DenseDoubleMatrix2D(A.rows(), A.columns());
        A_shrink.assign(weight).assign(DoubleFunctions.mult(-1.0D * epsilon)).assign(A, DoubleFunctions.plus).assign(DoubleFunctions.max(0.0D));
        A_shrinknegative.assign(weight).assign(DoubleFunctions.mult(epsilon)).assign(A, DoubleFunctions.plus).assign(DoubleFunctions.min(0.0D));
        A_shrink.assign(A_shrinknegative, DoubleFunctions.plus);
        return A_shrink;
    }

    public final DoubleMatrix2D stackToMatrix(ImageStack stack) {
        int rows = stack.getHeight();
        int columns = stack.getWidth();
        int nSlices = stack.getSize();
        DoubleMatrix2D stackMatrix = new DenseDoubleMatrix2D(rows * columns, nSlices);

        for(int k = 1; k <= nSlices; ++k) {
            for(int i = 0; i < stack.getWidth(); ++i) {
                for(int j = 0; j < stack.getHeight(); ++j) {
                    double temp1 = stack.getVoxel(i, j, k - 1);
                    stackMatrix.setQuick(i * rows + j, k - 1, temp1);
                }
            }
        }

        return stackMatrix;
    }

    public final ImageStack matrixToStack(DoubleMatrix2D matrix, int stackWidth, int stackHeight) {
        int stackSize = stackWidth * stackHeight;
        if (matrix.rows() != stackSize) {
            throw new EmptyStackException();
        } else {
            int stackZSlices = matrix.columns();
            new ImageStack(stackWidth, stackHeight, stackZSlices);
            ImagePlus imgPlus = NewImage.createFloatImage((String)null, stackWidth, stackHeight, stackZSlices, 2);
            ImageStack stack = imgPlus.getStack();

            for(int k = 1; k <= stackZSlices; ++k) {
                ImageProcessor ip = this.matrixToImage(matrix.viewColumn(k - 1).reshape(stackHeight, stackWidth));
                stack.setProcessor(ip, k);
            }

            return stack;
        }
    }

    public final ImageProcessor matrixToImage(DoubleMatrix2D matrix) {
        int imageWidth = matrix.columns();
        int imageHeight = matrix.rows();
        ImagePlus imgPlus = NewImage.createFloatImage((String)null, imageWidth, imageHeight, 1, 2);
        ImageProcessor ip = imgPlus.getProcessor();

        for(int i = 0; i < imageWidth; ++i) {
            for(int j = 0; j < imageHeight; ++j) {
                double temp1 = matrix.getQuick(i, j);
                ip.putPixelValue(j, i, temp1);
            }
        }

        return ip;
    }

    public final DoubleMatrix2D imageToMatrix(ImageProcessor ip) {
        int rows = ip.getHeight();
        int columns = ip.getWidth();
        DoubleMatrix2D matrix = new DenseDoubleMatrix2D(rows, columns);

        for(int i = 0; i < ip.getWidth(); ++i) {
            for(int j = 0; j < ip.getHeight(); ++j) {
                matrix.setQuick(j, i, (double)ip.getf(i, j));
            }
        }

        return matrix;
    }

    private class Baseline {
        public double[] basefluor;

        private Baseline() {
        }
    }

    private class DecomposedMatrix {
        public DoubleMatrix2D LowrankComponent;
        public DoubleMatrix2D SparseComponent;
        public DoubleMatrix1D Offset;
        public DoubleMatrix1D Coeff;

        private DecomposedMatrix() {
        }
    }

    private class Options {
        boolean lambda_auto;
        boolean lambdadark_auto;
        double lambda;
        double lambda_dark;
        boolean shadingEst;
        boolean darkfieldEst;
        boolean imageCorr;
        int driftOpt;

        private Options() {
            this.lambda_auto = true;
            this.lambdadark_auto = true;
            this.lambda = 0.5D;
            this.lambda_dark = 0.5D;
            this.shadingEst = true;
            this.darkfieldEst = false;
            this.imageCorr = false;
            this.driftOpt = 0;
        }
    }

    private static final class Parameters {
        private static final double epslon = 0.1D;
        private static final int reweightingIteration = 5;
        private static final int processingWidth = 128;
        private static final int processingHeight = 128;
        private static final double tolerance = 1.0E-6D;
        private static final double maxIter = 500.0D;

        private Parameters() {
        }
    }

    private class Shading {
        public ImagePlus flatfield;
        public ImagePlus darkfield;

        private Shading() {
        }
    }
}
