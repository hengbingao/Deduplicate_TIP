# Deduplicate_TIP
The methods for TIP-seq deduplication

![TIP-seq](https://github.com/hengbingao/Deduplicate_TIP/blob/main/png/Bulk_TIP-seq.png)


## **Install**

1. Clone the repository:

    ```bash
    git clone https://github.com/hengbingao/Deduplicate_TIP.git
    ```

2. Set the executable permissions:

    ```bash
    chmod +x $Deduplicate_TIP_PATH/bin/*
    chmod +x $Deduplicate_TIP_PATH/src/*
    ```

3. Add to environment:

    ```bash
    echo 'export PATH=$PATH:$Deduplicate_TIP_PATH/bin' >> ~/.bashrc
    source ~/.bashrc
    ```
## **Usage**

1. Uasge example:

    ```bash
    Deduplicate_TIP -i <input unsort.bam> -o <output deduplicated.bam> 
    ```

